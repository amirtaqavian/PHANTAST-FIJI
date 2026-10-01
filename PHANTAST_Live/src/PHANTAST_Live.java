/*
 * PHANTAST Live - Plugin for FIJI
 *
 * Re-implementation of PHANTAST (Jaccard et al., 2014) with:
 *  - sliders for sigma and epsilon
 *  - live confluency readout and coloured preview while adjusting
 *  - robust background preview (no crashes / no restart needed on extreme values)
 *  - confluency reported as a percentage
 *  - analysis area: whole image, user selection, or auto-detected circular field of view
 *  - image adjustments (brightness, contrast, clarity, gamma, CLAHE, median, smoothing,
 *    rolling-ball background, even illumination), shown live and optionally applied to detection
 *  - live yellow outline / green fill preview and live black & white mask window
 *  - tabbed settings window with saved presets (always opens on Default)
 *  - detection controls: min cell size, hole filling, grow/shrink, halo correction depth,
 *    exclusion of round bright cells, manual add/remove corrections
 * Uses only the core ImageJ API (no imglib2 dependency).
 *
 * Based on PHANTAST for FIJI:
 * Copyright (c) 2013, Nicolas Jaccard
 * Department of Biochemical Engineering, UCL
 * Centre for Mathematics and Physics in the Life Sciences and Experimental Biology, UCL
 * The British Heart Foundation
 * December 2017 modifications by Olivier Burri, BIOP, EPFL
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the distribution
 *     * Neither the names of the University College London or British Heart Foundation nor the names
 *       of its contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.Prefs;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.BackgroundSubtracter;
import ij.plugin.filter.EDM;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.RankFilters;
import ij.plugin.filter.ThresholdToSelection;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.TreeSet;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BooleanSupplier;

public class PHANTAST_Live implements PlugIn {

    // Slider ranges; typed values are clamped to these
    static final double SIGMA_MIN = 0.3, SIGMA_MAX = 10.0;
    static final double EPSILON_MIN = 0.001, EPSILON_MAX = 0.3;
    static final int MARGIN_MAX = 200;
    static final int HALO_FULL = 40;           // halo depth slider value meaning "no limit"
    static final int MIN_SIZE_MAX = 2000;
    static final int HOLE_MAX = 5000;
    static final int GROW_MAX = 10;

    // Analysis area modes
    static final int AREA_WHOLE = 0, AREA_SELECTION = 1, AREA_AUTO = 2;
    static final String[] AREA_NAMES = {"Whole image", "Selection drawn on image", "Auto-detect circular field"};

    // How detected cells are drawn in the live preview
    static final int STYLE_OUTLINE = 0, STYLE_FILL = 1, STYLE_BOTH = 2;
    static final String[] STYLE_NAMES = {"Yellow outline", "Green fill", "Yellow outline + green fill"};

    static final Color CELL_OUTLINE = Color.YELLOW;
    static final Color AREA_OUTLINE = Color.CYAN;
    static final byte FG = (byte) 255;
    static final String DEFAULT_PRESET = "Default";
    static final String PRESET_FILE = "PHANTAST_Live_presets.properties";

    // Kirsch compass kernels, projection cones and offsets from the original plugin
    static final float[][] KIRSCH = {
        {-3,-3,5,-3,0,5,-3,-3,5},
        {-3,5,5,-3,0,5,-3,-3,-3},
        {5,5,5,-3,0,-3,-3,-3,-3},
        {5,5,-3,5,0,-3,-3,-3,-3},
        {5,-3,-3,5,0,-3,5,-3,-3},
        {-3,-3,-3,5,0,-3,5,5,-3},
        {-3,-3,-3,-3,0,-3,5,5,5},
        {-3,-3,-3,-3,0,5,-3,5,5}};
    static final int[][] CONES = {{1,2,8},{2,1,3},{3,2,4},{4,3,5},{5,4,6},{6,5,7},{7,6,8},{8,1,7}};
    static final int[][] OFFSETS = {{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{-1,1},{0,1},{1,1}};

    // Detection settings (these determine confluency)
    private double sigma = 1.2;
    private Params params = Params.DEFAULT;
    private int areaMode = AREA_WHOLE;
    private int margin = 15;
    private Roi areaRoi;
    private int areaVersion;

    // Image adjustments: the sliders change the display instantly; they only reach
    // detection when "Apply to detection" is pressed (stored in applied)
    static final Display NEUTRAL = new Display(0, 1, 0, 1, 0, 0, 0, 0, false);
    private Display display = NEUTRAL;
    private Display applied = NEUTRAL;

    // Manual corrections, in the order they were made
    private final List<Edit> edits = new ArrayList<>();
    private int editVersion;

    // Preview settings
    private boolean preview = true;
    private int style = STYLE_OUTLINE;
    private boolean liveMask = false;
    // View only (not part of presets): hide the detection drawing while the Image adjustments tab is open
    private volatile boolean hideDetection = true;
    private volatile boolean onAdjustTab;
    static final int ADJUST_TAB = 1;
    private JCheckBox bHideDetection;

    // Outputs when OK is clicked
    private boolean outputTable = true;
    private boolean outputOutline = true;
    private boolean outputMask = true;
    private boolean allSlices = false;

    private ImagePlus imp;
    private Overlay originalOverlay;
    private Roi originalRoi;
    private JDialog dialog;
    private JLabel status;
    private JLabel appliedLabel;
    private JLabel editsLabel;
    private JLabel areaLabel;
    private JLabel imageLabel;
    private JButton setAreaButton;
    private String lastDetectionKey = "";
    private String lastDisplayKey = "";

    private final AtomicInteger detectGen = new AtomicInteger();
    private final AtomicInteger displayGen = new AtomicInteger();
    private final ExecutorService detectWorker = daemonThread("PHANTAST Live detection");
    private final ExecutorService displayWorker = daemonThread("PHANTAST Live display");

    // Latest preview layers, written by the workers and drawn on the event thread
    private volatile ImageRoi layerView;
    private volatile ImageRoi layerFill;
    private volatile Roi layerOutline;
    private volatile Roi layerArea;
    private volatile ByteProcessor layerMask;
    private ImagePlus liveMaskImp;

    // Detection cache, only touched by the detection worker
    private ImagePlus dImp, vImp;
    private int dSlice = -1;
    private FloatProcessor dImage;
    private String dAdjustKey;
    private FloatProcessor dInput;
    private double dSigma = Double.NaN;
    private float[] dContrast;
    private byte[] dDirections;
    private int dAutoMargin = -1;
    private Area dAutoArea;

    // Display cache, only touched by the display worker
    private int vSlice = -1;
    private FloatProcessor vRaw;
    private String vBaseKey;
    private float[] vBase;
    private float[] vBaseBlur;
    private double vMean;

    // Window controls
    private SliderField fMargin, fBrightness, fContrast, fClarity, fGamma, fClahe, fMedian, fSmooth, fBackground;
    private SliderField fSigma, fEpsilon, fHalo, fMinSize, fHoles, fGrow;
    private JComboBox<String> cArea, cStyle, cPreset;
    private JCheckBox bFlatten, bFillAll, bRound, bPreview, bLiveMask, bTable, bOutline, bMask, bAllSlices;
    private JLabel presetState;

    // Batch tab
    private JTextField batchFolder;
    private JCheckBox bSubfolders, bSaveMasks, bSaveOutlines;
    private JComboBox<String> cBatchPreset;
    private JButton runBatchButton, stopBatchButton, openResultsButton;
    private javax.swing.JProgressBar batchProgress;
    private JLabel batchStatus;
    private volatile boolean batchCancel;
    private Thread batchThread;
    private File lastBatchOutput;
    private boolean updatingControls;
    private String loadedPresetName = DEFAULT_PRESET;
    private String loadedPresetValues;
    private final CountDownLatch closed = new CountDownLatch(1);
    private volatile boolean okPressed;

    /** Thrown to abandon a preview that has been superseded by newer settings. */
    static class Cancelled extends RuntimeException {
        Cancelled() { super(null, null, false, false); }
    }

    /** Analysis area as a pixel mask (null = whole image) plus an outline to draw. */
    static final class Area {
        final byte[] mask;
        final Roi outline;
        final String error;
        Area(byte[] mask, Roi outline, String error) {
            this.mask = mask;
            this.outline = outline;
            this.error = error;
        }
    }

    /** One manual correction: a drawn region added to or removed from the cells. */
    static final class Edit {
        final int slice;
        final Roi roi;
        final boolean add;
        Edit(int slice, Roi roi, boolean add) {
            this.slice = slice;
            this.roi = roi;
            this.add = add;
        }
    }

    /** Detection settings applied after the local-contrast step. */
    static final class Params {
        static final Params DEFAULT = new Params(0.03, HALO_FULL, 100, 25, false, 0, false);
        final double epsilon;
        final int haloDepth, minSize, maxHole, grow;
        final boolean fillAll, excludeRound;

        Params(double epsilon, int haloDepth, int minSize, int maxHole, boolean fillAll, int grow,
               boolean excludeRound) {
            this.epsilon = clamp(epsilon, EPSILON_MIN, EPSILON_MAX);
            this.haloDepth = (int) clamp(haloDepth, 0, HALO_FULL);
            this.minSize = (int) clamp(minSize, 0, MIN_SIZE_MAX);
            this.maxHole = (int) clamp(maxHole, 0, HOLE_MAX);
            this.fillAll = fillAll;
            this.grow = (int) clamp(grow, -GROW_MAX, GROW_MAX);
            this.excludeRound = excludeRound;
        }

        String key() {
            return epsilon + "|" + haloDepth + "|" + minSize + "|" + maxHole + "|" + fillAll + "|" + grow + "|" + excludeRound;
        }
    }

    /**
     * Image adjustments. Shown live on screen; used for detection only once applied.
     * The "base" steps need filtering and are cached; the rest is a fast per-pixel pass.
     */
    static final class Display {
        static final double BRIGHTNESS_MAX = 0.5, CONTRAST_MIN = 0.2, CONTRAST_MAX = 3.0;
        static final double CLARITY_MAX = 10.0, GAMMA_MIN = 0.2, GAMMA_MAX = 5.0, CLAHE_MAX = 10.0;
        static final double MEDIAN_MAX = 10.0, DENOISE_MAX = 5.0, BACKGROUND_MAX = 300.0;
        static final double CLARITY_RADIUS = 8.0;
        static final int CLAHE_TILE = 128;
        final double brightness, contrast, clarity, gamma, clahe, median, denoise, background;
        final boolean flatten;

        Display(double brightness, double contrast, double clarity, double gamma, double clahe, double median,
                double denoise, double background, boolean flatten) {
            this.brightness = clamp(brightness, -BRIGHTNESS_MAX, BRIGHTNESS_MAX);
            this.contrast = clamp(contrast, CONTRAST_MIN, CONTRAST_MAX);
            this.clarity = clamp(clarity, 0, CLARITY_MAX);
            this.gamma = clamp(gamma, GAMMA_MIN, GAMMA_MAX);
            this.clahe = clamp(clahe, 0, CLAHE_MAX);
            this.median = clamp(median, 0, MEDIAN_MAX);
            this.denoise = clamp(denoise, 0, DENOISE_MAX);
            this.background = clamp(background, 0, BACKGROUND_MAX);
            this.flatten = flatten;
        }

        boolean neutral() {
            return key().equals(NEUTRAL_KEY);
        }

        String key() {
            return baseKey() + "|" + brightness + "|" + contrast + "|" + clarity + "|" + gamma;
        }

        String baseKey() {
            return background + "|" + flatten + "|" + median + "|" + denoise + "|" + clahe;
        }

        String describe() {
            if (neutral()) return "None";
            List<String> parts = new ArrayList<>();
            if (brightness != 0) parts.add(String.format("brightness %+.0f%%", brightness * 100));
            if (contrast != 1) parts.add(String.format("contrast %.2f", contrast));
            if (clarity != 0) parts.add(String.format("clarity %.2f", clarity));
            if (gamma != 1) parts.add(String.format("gamma %.2f", gamma));
            if (clahe != 0) parts.add(String.format("CLAHE %.1f", clahe));
            if (median != 0) parts.add(String.format("median %.1f px", median));
            if (denoise != 0) parts.add(String.format("smooth %.1f px", denoise));
            if (background != 0) parts.add(String.format("background %.0f px", background));
            if (flatten) parts.add("even illumination");
            return String.join(", ", parts);
        }

        /** Filtering steps (cached between slider moves). Input and output are 0..1. */
        FloatProcessor buildBase(FloatProcessor view) {
            FloatProcessor base = view;
            if (background > 0) base = rollingBallFlatten(base, background);
            if (flatten) base = flatten(base);
            if (median > 0) {
                base = (FloatProcessor) base.duplicate();
                new RankFilters().rank(base, median, RankFilters.MEDIAN);
            }
            if (denoise > 0) {
                base = (FloatProcessor) base.duplicate();
                new GaussianBlur().blurGaussian(base, denoise, denoise, 0.002);
            }
            if (clahe > 0) base = clahe(base, clahe, CLAHE_TILE);
            return base;
        }

        /** The adjusted image exactly as shown on screen (0..1), used when applied to detection. */
        FloatProcessor apply(FloatProcessor view) {
            if (neutral()) return view;
            FloatProcessor base = buildBase(view);
            float[] b = (float[]) base.getPixels();
            float[] bb = null;
            if (clarity > 0) {
                FloatProcessor blur = new FloatProcessor(base.getWidth(), base.getHeight(), b.clone());
                new GaussianBlur().blurGaussian(blur, CLARITY_RADIUS, CLARITY_RADIUS, 0.01);
                bb = (float[]) blur.getPixels();
            }
            byte[] out = renderPixels(b, bb, mean(b), this);
            float[] f = new float[out.length];
            for (int i = 0; i < out.length; i++) f[i] = (out[i] & 0xff) / 255f;
            return new FloatProcessor(base.getWidth(), base.getHeight(), f);
        }
    }

    static final String NEUTRAL_KEY = new Display(0, 1, 0, 1, 0, 0, 0, 0, false).key();

    /** Result of one detection run. */
    static final class Detection {
        final ByteProcessor mask;
        final double confluency;
        final Area area;
        Detection(ByteProcessor mask, double confluency, Area area) {
            this.mask = mask;
            this.confluency = confluency;
            this.area = area;
        }
    }

    private static ExecutorService daemonThread(String name) {
        return Executors.newSingleThreadExecutor(r -> {
            Thread t = new Thread(r, name);
            t.setDaemon(true);
            return t;
        });
    }

    // ---------------------------------------------------------------- run

    @Override
    public void run(String arg) {
        // Macro / batch use: run directly with a saved preset, e.g. run("PHANTAST Live", "preset=[My cells]")
        String options = Macro.getOptions();
        if (options != null) {
            imp = IJ.getImage(); // a macro needs an open image
            if (imp == null) return;
            originalOverlay = imp.getOverlay();
            originalRoi = imp.getRoi();
            if (originalRoi != null && originalRoi.isArea()) areaRoi = (Roi) originalRoi.clone();
            String name = Macro.getValue(options, "preset", DEFAULT_PRESET).trim();
            if (!DEFAULT_PRESET.equalsIgnoreCase(name)) {
                String values = loadPresets().get(name);
                if (values == null) {
                    IJ.error("PHANTAST Live", "Preset not found: " + name);
                    return;
                }
                fieldsFromMap(parse(values));
            }
            finalRun();
            return;
        }

        // The window opens with or without an image; an image can be opened later
        final ImagePlus start = activeImage();
        if (start != null && start.getRoi() != null && start.getRoi().isArea()) areaMode = AREA_SELECTION;
        try {
            SwingUtilities.invokeAndWait(() -> {
                buildDialog();
                attachImage(start);
                placeDialog();
                dialog.setVisible(true);
            });
        } catch (Exception e) {
            IJ.handleException(e);
            return;
        }
        ImagePlus.addImageListener(imageListener);
        try {
            closed.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        ImagePlus.removeImageListener(imageListener);

        detectGen.incrementAndGet(); // abandon any running preview
        displayGen.incrementAndGet();
        detectWorker.shutdown();
        displayWorker.shutdown();
        try {
            SwingUtilities.invokeAndWait(() -> {
                closeLiveMask();
                restoreImage();
            });
        } catch (Exception e) {
            IJ.handleException(e);
        }
        if (!okPressed || imp == null) return;
        try {
            finalRun();
        } catch (Throwable t) {
            IJ.handleException(t);
        }
    }

    // ---------------------------------------------------------------- which image the window works on

    /** Follows images being opened and closed while the window is open. */
    private final ij.ImageListener imageListener = new ij.ImageListener() {
        @Override
        public void imageOpened(ImagePlus im) {
            SwingUtilities.invokeLater(() -> {
                if (dialog != null && imp == null && !isOwnWindow(im)) attachImage(im);
            });
        }

        @Override
        public void imageClosed(ImagePlus im) {
            SwingUtilities.invokeLater(() -> {
                if (dialog == null || im != imp) return;
                imp = null; // closed: nothing to restore
                ImagePlus next = activeImage();
                attachImage(next == im ? null : next);
            });
        }

        @Override
        public void imageUpdated(ImagePlus im) {
        }
    };

    /** The image in the front window, ignoring this plugin's own mask window. */
    private ImagePlus activeImage() {
        ImagePlus im = ij.WindowManager.getCurrentImage();
        return im == null || isOwnWindow(im) ? null : im;
    }

    private boolean isOwnWindow(ImagePlus im) {
        return im == liveMaskImp || (im.getTitle() != null && im.getTitle().endsWith(" - live mask"));
    }

    /** Puts the image back as it was (overlay and selection). */
    private void restoreImage() {
        if (imp == null) return;
        imp.setOverlay(originalOverlay);
        if (originalRoi != null) imp.setRoi(originalRoi);
    }

    /** Makes the window work on {@code im} (or on no image); drawn area and manual edits belong to one image. */
    private void attachImage(ImagePlus im) {
        if (im != null && im == imp) return;
        restoreImage();
        closeLiveMask();
        imp = im;
        edits.clear();
        editVersion++;
        areaRoi = null;
        areaVersion++;
        layerView = null; layerFill = null; layerOutline = null; layerArea = null; layerMask = null;
        originalOverlay = null;
        originalRoi = null;
        if (imp != null) {
            originalOverlay = imp.getOverlay();
            originalRoi = imp.getRoi();
            if (originalRoi != null && originalRoi.isArea()) {
                areaRoi = (Roi) originalRoi.clone();
                imp.deleteRoi(); // free the image selection for manual edits
            }
        }
        if (imageLabel != null) {
            imageLabel.setText(imp == null ? "Image: none open - open an image (File > Open) to see the live preview"
                    : "Image: " + imp.getTitle() + (imp.getStackSize() > 1 ? "  (" + imp.getStackSize() + " slices)" : ""));
            bAllSlices.setEnabled(imp != null && imp.getStackSize() > 1);
            updateAreaLabel();
            updateEditsLabel();
        }
        lastDetectionKey = detectionKey();
        lastDisplayKey = display.key();
        scheduleDetection();
        scheduleDisplay();
    }

    private void useActiveImage() {
        ImagePlus im = activeImage();
        if (im == null) {
            setStatus("No image open - open one with File > Open");
        } else if (im == imp) {
            setStatus("Already using " + im.getTitle());
        } else {
            attachImage(im);
        }
    }

    // ---------------------------------------------------------------- batch tab

    static final String BATCH_CURRENT = "Current settings (this window)";

    private void refreshBatchPresets() {
        if (cBatchPreset == null) return;
        Object previous = cBatchPreset.getSelectedItem();
        cBatchPreset.removeAllItems();
        cBatchPreset.addItem(BATCH_CURRENT);
        cBatchPreset.addItem(DEFAULT_PRESET);
        for (String name : new TreeSet<>(loadPresets().keySet())) cBatchPreset.addItem(name);
        if (previous != null) cBatchPreset.setSelectedItem(previous);
    }

    private void chooseBatchFolder() {
        javax.swing.JFileChooser fc = new javax.swing.JFileChooser(batchFolder.getText().trim());
        fc.setFileSelectionMode(javax.swing.JFileChooser.DIRECTORIES_ONLY);
        fc.setDialogTitle("Choose the folder with your images");
        if (fc.showOpenDialog(dialog) == javax.swing.JFileChooser.APPROVE_OPTION)
            batchFolder.setText(fc.getSelectedFile().getAbsolutePath());
    }

    private void startBatch() {
        if (batchThread != null) return;
        File dir = new File(batchFolder.getText().trim());
        if (batchFolder.getText().trim().isEmpty() || !dir.isDirectory()) {
            batchStatus.setText("Choose a folder first");
            return;
        }
        String choice = String.valueOf(cBatchPreset.getSelectedItem());
        Map<String, String> settings;
        if (BATCH_CURRENT.equals(choice)) {
            settings = currentSettings();
        } else if (DEFAULT_PRESET.equals(choice)) {
            settings = defaultSettings();
        } else {
            String stored = loadPresets().get(choice);
            if (stored == null) {
                batchStatus.setText("Preset not found: " + choice);
                return;
            }
            settings = parse(stored);
        }
        final Roi area = areaRoi == null ? null : (Roi) areaRoi.clone();
        final boolean recurse = bSubfolders.isSelected(), masks = bSaveMasks.isSelected(),
                outlines = bSaveOutlines.isSelected();
        batchCancel = false;
        runBatchButton.setEnabled(false);
        stopBatchButton.setEnabled(true);
        openResultsButton.setEnabled(false);
        batchProgress.setValue(0);
        batchProgress.setString("Starting...");
        batchStatus.setText("Settings: " + choice);
        batchThread = new Thread(() -> {
            BatchResult r;
            try {
                r = runBatch(dir, recurse, settings, area, masks, outlines, () -> batchCancel,
                        (done, total, message) -> SwingUtilities.invokeLater(() -> {
                            batchProgress.setMaximum(total);
                            batchProgress.setValue(done);
                            batchProgress.setString(done + " / " + total + " images");
                            batchStatus.setText(message);
                        }));
            } catch (Throwable t) {
                r = new BatchResult();
                r.error = "Batch failed: " + t;
                IJ.log("PHANTAST Live batch error: " + t);
            }
            final BatchResult result = r;
            SwingUtilities.invokeLater(() -> batchFinished(result));
        }, "PHANTAST Live batch");
        batchThread.setDaemon(true);
        batchThread.start();
    }

    private void batchFinished(BatchResult r) {
        batchThread = null;
        if (runBatchButton == null) return;
        runBatchButton.setEnabled(true);
        stopBatchButton.setEnabled(false);
        if (r.error != null) {
            batchProgress.setString("");
            batchStatus.setText("<html>" + r.error + "</html>");
            return;
        }
        lastBatchOutput = r.outputDir;
        openResultsButton.setEnabled(true);
        batchProgress.setString((r.cancelled ? "Stopped: " : "Done: ") + r.images + " image(s)"
                + (r.failed > 0 ? ", " + r.failed + " could not be read" : ""));
        batchStatus.setText("Saved in " + r.outputDir.getName());
        if (r.table != null && r.table.size() > 0) r.table.show("PHANTAST Batch Results");
    }

    private void openBatchResults() {
        if (lastBatchOutput == null) return;
        try {
            java.awt.Desktop.getDesktop().open(lastBatchOutput);
        } catch (Exception e) {
            batchStatus.setText("Results are in " + lastBatchOutput.getAbsolutePath());
        }
    }

    private void finish(boolean ok) {
        if (dialog == null) return;
        if (ok && imp == null) {
            setStatus("Open an image to measure it (or use Cancel to close)");
            return;
        }
        batchCancel = true; // a running batch stops after the current image
        okPressed = ok;
        dialog.dispose();
        dialog = null;
        closed.countDown();
    }

    // ---------------------------------------------------------------- window

    /** A labelled slider with a text box; typed values are clamped to the range. */
    final class SliderField {
        final double min, max, step;
        final int decimals;
        final JLabel label;
        final JSlider slider;
        final JTextField text;
        double value;
        private boolean syncing;

        SliderField(String name, double min, double max, double value, double step) {
            this.min = min;
            this.max = max;
            this.step = step;
            this.decimals = step >= 1 ? 0 : step >= 0.1 ? 1 : step >= 0.01 ? 2 : 3;
            label = new JLabel(name);
            slider = new JSlider(0, (int) Math.round((max - min) / step), 0);
            slider.setPreferredSize(new Dimension(190, slider.getPreferredSize().height));
            text = new JTextField(6);
            text.setHorizontalAlignment(JTextField.RIGHT);
            setValue(value);
            slider.addChangeListener(e -> {
                if (syncing) return;
                setValue(min + slider.getValue() * step);
                controlsChanged();
            });
            text.addActionListener(e -> commitText());
            text.addFocusListener(new FocusAdapter() {
                @Override
                public void focusLost(FocusEvent e) {
                    commitText();
                }
            });
        }

        void commitText() {
            try {
                double v = Double.parseDouble(text.getText().trim());
                if (Double.isNaN(v) || Double.isInfinite(v)) throw new NumberFormatException();
                double before = value;
                setValue(v);
                if (value != before) controlsChanged();
            } catch (NumberFormatException ex) {
                setValue(value); // put the last valid value back
            }
        }

        void setValue(double v) {
            double f = Math.pow(10, decimals);
            value = Math.round(clamp(v, min, max) * f) / f;
            syncing = true;
            slider.setValue((int) Math.round((value - min) / step));
            text.setText(IJ.d2s(value, decimals));
            syncing = false;
        }
    }

    /** Simple label / slider / value form. */
    static final class Form {
        final JPanel panel = new JPanel(new GridBagLayout());
        int row;

        Form() {
            panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        }

        void slider(SliderField f) {
            GridBagConstraints c = new GridBagConstraints();
            c.gridy = row++;
            c.insets = new Insets(3, 0, 3, 8);
            c.anchor = GridBagConstraints.WEST;
            c.gridx = 0;
            panel.add(f.label, c);
            c.gridx = 1;
            c.weightx = 1;
            c.fill = GridBagConstraints.HORIZONTAL;
            panel.add(f.slider, c);
            c.gridx = 2;
            c.weightx = 0;
            c.fill = GridBagConstraints.NONE;
            c.insets = new Insets(3, 0, 3, 0);
            panel.add(f.text, c);
        }

        void labelled(String name, Component comp) {
            GridBagConstraints c = new GridBagConstraints();
            c.gridy = row++;
            c.insets = new Insets(3, 0, 3, 8);
            c.anchor = GridBagConstraints.WEST;
            c.gridx = 0;
            panel.add(new JLabel(name), c);
            c.gridx = 1;
            c.gridwidth = 2;
            panel.add(comp, c);
        }

        void full(Component comp) {
            GridBagConstraints c = new GridBagConstraints();
            c.gridy = row++;
            c.gridx = 0;
            c.gridwidth = 3;
            c.insets = new Insets(4, 0, 4, 0);
            c.anchor = GridBagConstraints.WEST;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.weightx = 1;
            panel.add(comp, c);
        }

        JComponent done() {
            JPanel wrap = new JPanel(new BorderLayout());
            wrap.add(panel, BorderLayout.NORTH);
            return wrap;
        }
    }

    private static JLabel help(String html) {
        JLabel l = new JLabel("<html><div style='width:380px'>" + html + "</div></html>");
        l.setFont(l.getFont().deriveFont(Font.ITALIC, 11f));
        l.setForeground(new Color(90, 90, 90));
        return l;
    }

    private static JButton button(String label, Runnable action) {
        JButton b = new JButton(label);
        b.addActionListener(e -> action.run());
        return b;
    }

    private JCheckBox checkbox(String label, boolean value) {
        JCheckBox b = new JCheckBox(label, value);
        b.addItemListener(e -> controlsChanged());
        return b;
    }

    private JComboBox<String> combo(String[] items, int selected) {
        JComboBox<String> c = new JComboBox<>(items);
        c.setSelectedIndex(selected);
        c.addActionListener(e -> controlsChanged());
        return c;
    }

    /** Builds and shows the settings window, then starts the live preview; runs on the event thread. */
    void buildDialog() {
        createDialog();
    }

    /** Builds the settings window without showing it. */
    void createDialog() {
        loadedPresetValues = serialize(currentSettings());
        Display d = display;

        // 1. Area
        Form area = new Form();
        cArea = combo(AREA_NAMES, areaMode);
        area.labelled("Analysis area", cArea);
        fMargin = new SliderField("Edge margin (px, auto area)", 0, MARGIN_MAX, margin, 1);
        area.slider(fMargin);
        setAreaButton = button("Set area from selection", this::setAreaFromSelection);
        area.full(setAreaButton);
        areaLabel = new JLabel(" ");
        area.full(areaLabel);
        area.full(help("<b>Whole image:</b> analyse everything.<br><b>Selection:</b> draw the area on the image "
                + "(e.g. oval tool) and press <i>Set area from selection</i>.<br><b>Auto-detect:</b> finds the bright "
                + "circle of photos taken through the eyepiece; the margin trims its edge."));

        // 2. Image adjustments
        Form adjust = new Form();
        fBrightness = new SliderField("Brightness (%)", -50, 50, d.brightness * 100, 1);
        fContrast = new SliderField("Contrast (1 = unchanged)", Display.CONTRAST_MIN, Display.CONTRAST_MAX, d.contrast, 0.05);
        fClarity = new SliderField("Clarity (0 = off)", 0, Display.CLARITY_MAX, d.clarity, 0.05);
        fGamma = new SliderField("Gamma (<1 brightens dark areas)", Display.GAMMA_MIN, Display.GAMMA_MAX, d.gamma, 0.05);
        fClahe = new SliderField("Local contrast / CLAHE (0 = off)", 0, Display.CLAHE_MAX, d.clahe, 0.1);
        fMedian = new SliderField("Median smoothing (px, 0 = off)", 0, Display.MEDIAN_MAX, d.median, 0.5);
        fSmooth = new SliderField("Gaussian smoothing (px, 0 = off)", 0, Display.DENOISE_MAX, d.denoise, 0.1);
        fBackground = new SliderField("Background correction (px, 0 = off)", 0, Display.BACKGROUND_MAX, d.background, 5);
        for (SliderField f : new SliderField[] {fBrightness, fContrast, fClarity, fGamma, fClahe, fMedian, fSmooth, fBackground})
            adjust.slider(f);
        bFlatten = checkbox("Even out uneven illumination", d.flatten);
        adjust.full(bFlatten);
        bHideDetection = new JCheckBox("Hide cell detection while on this tab (see the cells clearly)", hideDetection);
        bHideDetection.addItemListener(e -> {
            hideDetection = bHideDetection.isSelected();
            rebuildOverlay();
        });
        adjust.full(bHideDetection);
        JPanel adjustButtons = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        adjustButtons.add(button("Apply to detection", this::applyAdjustments));
        adjustButtons.add(Box.createHorizontalStrut(8));
        adjustButtons.add(button("Default", this::resetAdjustments));
        adjust.full(adjustButtons);
        adjust.full(help("Changes show on the image instantly. Press <i>Apply to detection</i> to detect cells on the "
                + "adjusted image; <i>Default</i> resets these sliders and detection to the original image."));

        // 3. Cell detection (preview display options on top)
        Form detect = new Form();
        bPreview = checkbox("Live preview", preview);
        cStyle = combo(STYLE_NAMES, style);
        bLiveMask = checkbox("Live black & white mask window", liveMask);
        JPanel previewRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        previewRow.add(bPreview);
        previewRow.add(Box.createHorizontalStrut(12));
        previewRow.add(new JLabel("Show cells as "));
        previewRow.add(cStyle);
        detect.full(previewRow);
        detect.full(bLiveMask);
        detect.full(new javax.swing.JSeparator());
        fSigma = new SliderField("Sigma", SIGMA_MIN, SIGMA_MAX, sigma, 0.1);
        fEpsilon = new SliderField("Epsilon", EPSILON_MIN, EPSILON_MAX, params.epsilon, 0.001);
        fHalo = new SliderField("Halo correction (px, 0=off, 40=full)", 0, HALO_FULL, params.haloDepth, 1);
        fMinSize = new SliderField("Minimum cell size (px)", 0, MIN_SIZE_MAX, params.minSize, 10);
        fHoles = new SliderField("Fill holes up to (px)", 0, HOLE_MAX, params.maxHole, 5);
        fGrow = new SliderField("Grow (+) / shrink (-) cells (px)", -GROW_MAX, GROW_MAX, params.grow, 1);
        for (SliderField f : new SliderField[] {fSigma, fEpsilon, fHalo, fMinSize, fHoles, fGrow}) detect.slider(f);
        bFillAll = checkbox("Fill all holes inside cells", params.fillAll);
        bRound = checkbox("Exclude round bright cells", params.excludeRound);
        detect.full(bFillAll);
        detect.full(bRound);

        // 4. Manual correction
        Form manual = new Form();
        manual.full(help("Draw a selection on the image (freehand, oval, polygon...), then press a button. "
                + "Edits are applied after automatic detection."));
        JPanel editButtons = new JPanel(new GridLayout(2, 2, 6, 6));
        editButtons.add(button("Add selection to cells", () -> addEdit(true)));
        editButtons.add(button("Remove selection", () -> addEdit(false)));
        editButtons.add(button("Undo last edit", this::undoEdit));
        editButtons.add(button("Clear all edits", this::clearEdits));
        manual.full(editButtons);
        editsLabel = new JLabel("Manual edits: none");
        manual.full(editsLabel);

        // 5. Output
        Form output = new Form();
        JLabel outHeader = new JLabel("When you click OK:");
        outHeader.setFont(outHeader.getFont().deriveFont(Font.BOLD));
        output.full(outHeader);
        bTable = checkbox("Add confluency to Results table", outputTable);
        bOutline = checkbox("Draw yellow cell outline on image", outputOutline);
        bMask = checkbox("Create black & white mask image", outputMask);
        bAllSlices = checkbox("Process all slices of the stack", allSlices);
        bAllSlices.setEnabled(imp != null && imp.getStackSize() > 1);
        for (JCheckBox b : new JCheckBox[] {bTable, bOutline, bMask, bAllSlices}) output.full(b);

        // 6. Batch process (not part of presets)
        Form batch = new Form();
        batchFolder = new JTextField(22);
        JPanel folderRow = new JPanel(new BorderLayout(6, 0));
        folderRow.add(batchFolder, BorderLayout.CENTER);
        folderRow.add(button("Choose...", this::chooseBatchFolder), BorderLayout.EAST);
        batch.labelled("Folder", folderRow);
        bSubfolders = new JCheckBox("Include subfolders");
        batch.full(bSubfolders);
        cBatchPreset = new JComboBox<>();
        refreshBatchPresets();
        batch.labelled("Settings", cBatchPreset);
        bSaveMasks = new JCheckBox("Save a black & white mask for each image", true);
        bSaveOutlines = new JCheckBox("Save each image with the yellow cell outline", true);
        batch.full(bSaveMasks);
        batch.full(bSaveOutlines);
        JPanel batchButtons = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        runBatchButton = button("Run batch", this::startBatch);
        stopBatchButton = button("Stop", () -> batchCancel = true);
        stopBatchButton.setEnabled(false);
        openResultsButton = button("Open results folder", this::openBatchResults);
        openResultsButton.setEnabled(false);
        batchButtons.add(runBatchButton);
        batchButtons.add(Box.createHorizontalStrut(6));
        batchButtons.add(stopBatchButton);
        batchButtons.add(Box.createHorizontalStrut(6));
        batchButtons.add(openResultsButton);
        batch.full(batchButtons);
        batchProgress = new javax.swing.JProgressBar(0, 1);
        batchProgress.setStringPainted(true);
        batchProgress.setString("");
        batch.full(batchProgress);
        batchStatus = new JLabel(" ");
        batch.full(batchStatus);
        batch.full(help("Every image in the folder is analysed with the chosen settings (manual corrections are not "
                + "used). Results are saved in a new folder <i>PHANTAST_Live_results_&lt;date&gt;</i> inside it: "
                + "<i>confluency_results.csv</i>, <i>settings_used.txt</i>, and the masks / outline images."));

        JTabbedPane tabs = new JTabbedPane();
        tabs.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT); // one row; never reshuffles on click
        tabs.addTab("1. Area", area.done());
        tabs.addTab("2. Image adjustments", adjust.done());
        tabs.addTab("3. Cell detection", detect.done());
        tabs.addTab("4. Manual correction", manual.done());
        tabs.addTab("5. Output", output.done());
        tabs.addTab("6. Batch process", batch.done());
        tabs.addChangeListener(e -> {
            onAdjustTab = tabs.getSelectedIndex() == ADJUST_TAB;
            rebuildOverlay();
        });

        // Preset bar
        JPanel presetBar = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        presetBar.add(new JLabel("Preset:"));
        cPreset = new JComboBox<>();
        cPreset.setPrototypeDisplayValue("XXXXXXXXXXXXXXXXXXXXXX");
        refreshPresetList(DEFAULT_PRESET);
        cPreset.addActionListener(e -> {
            if (updatingControls) return;
            Object sel = cPreset.getSelectedItem();
            if (sel != null) loadPreset(sel.toString());
        });
        presetBar.add(cPreset);
        presetBar.add(button("Save preset", this::savePreset));
        presetBar.add(button("Delete", this::deletePreset));
        presetState = new JLabel(" ");
        presetState.setForeground(new Color(180, 90, 0));
        presetBar.add(presetState);

        // Status and buttons
        status = new JLabel("Confluency: calculating...");
        status.setFont(status.getFont().deriveFont(Font.BOLD, 16f));
        appliedLabel = new JLabel(" ");
        appliedLabel.setFont(appliedLabel.getFont().deriveFont(Font.ITALIC, 11f));
        JPanel statusPanel = new JPanel();
        statusPanel.setLayout(new BoxLayout(statusPanel, BoxLayout.Y_AXIS));
        statusPanel.add(status);
        statusPanel.add(appliedLabel);
        JPanel okCancel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 6, 0));
        JButton ok = button("OK", () -> finish(true));
        okCancel.add(button("Cancel", () -> finish(false)));
        okCancel.add(ok);
        JPanel bottom = new JPanel(new BorderLayout(10, 0));
        bottom.setBorder(BorderFactory.createEmptyBorder(8, 2, 0, 2));
        bottom.add(statusPanel, BorderLayout.CENTER);
        bottom.add(okCancel, BorderLayout.EAST);

        JPanel root = new JPanel(new BorderLayout(0, 8));
        root.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
        imageLabel = new JLabel("Image: none open");
        JPanel imageBar = new JPanel(new BorderLayout(8, 0));
        imageBar.add(imageLabel, BorderLayout.CENTER);
        imageBar.add(button("Use active image", this::useActiveImage), BorderLayout.EAST);
        JPanel topBars = new JPanel(new BorderLayout(0, 6));
        topBars.add(imageBar, BorderLayout.NORTH);
        topBars.add(presetBar, BorderLayout.SOUTH);
        root.add(topBars, BorderLayout.NORTH);
        root.add(tabs, BorderLayout.CENTER);
        root.add(bottom, BorderLayout.SOUTH);

        dialog = new JDialog(IJ.getInstance(), "PHANTAST Live", false);
        dialog.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
        dialog.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                finish(false);
            }
        });
        dialog.setContentPane(root);
        dialog.getRootPane().setDefaultButton(ok);
        dialog.pack();
        widenForTabs(tabs);
        placeDialog();
        updateAreaLabel();
        updateEditsLabel();
        updateAppliedLabel();
        updatePresetState();
    }

    /** Widens the window until every tab fits side by side, and keeps that as the minimum size. */
    private void widenForTabs(JTabbedPane tabs) {
        int last = tabs.getTabCount() - 1;
        for (int i = 0; i < 60; i++) {
            dialog.validate();
            Rectangle r = tabs.getUI().getTabBounds(tabs, last);
            if (r != null && r.x + r.width + 8 <= tabs.getWidth()) break;
            dialog.setSize(dialog.getWidth() + 20, dialog.getHeight());
        }
        dialog.setMinimumSize(dialog.getSize());
    }

    /** Puts the window next to the image if there is room, otherwise keeps it on screen. */
    private void placeDialog() {
        Rectangle screen = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
        Dimension size = dialog.getSize();
        int x = screen.x + (screen.width - size.width) / 2, y = screen.y + (screen.height - size.height) / 2;
        if (imp != null && imp.getWindow() != null) {
            Rectangle win = imp.getWindow().getBounds();
            x = win.x + win.width + 8;
            y = win.y;
            if (x + size.width > screen.x + screen.width) x = screen.x + screen.width - size.width;
        }
        x = Math.max(screen.x, x);
        y = Math.max(screen.y, Math.min(y, screen.y + screen.height - size.height));
        dialog.setLocation(x, y);
    }

    /** Reads every control into the settings fields. */
    private void readControls() {
        sigma = fSigma.value;
        params = new Params(fEpsilon.value, (int) fHalo.value, (int) fMinSize.value, (int) fHoles.value,
                bFillAll.isSelected(), (int) fGrow.value, bRound.isSelected());
        areaMode = cArea.getSelectedIndex();
        margin = (int) fMargin.value;
        display = new Display(fBrightness.value / 100.0, fContrast.value, fClarity.value, fGamma.value,
                fClahe.value, fMedian.value, fSmooth.value, fBackground.value, bFlatten.isSelected());
        preview = bPreview.isSelected();
        style = cStyle.getSelectedIndex();
        liveMask = bLiveMask.isSelected();
        outputTable = bTable.isSelected();
        outputOutline = bOutline.isSelected();
        outputMask = bMask.isSelected();
        allSlices = bAllSlices.isSelected();
    }

    /** Puts the settings fields into the controls without triggering updates. */
    private void syncControls() {
        updatingControls = true;
        try {
            fSigma.setValue(sigma);
            fEpsilon.setValue(params.epsilon);
            fHalo.setValue(params.haloDepth);
            fMinSize.setValue(params.minSize);
            fHoles.setValue(params.maxHole);
            fGrow.setValue(params.grow);
            bFillAll.setSelected(params.fillAll);
            bRound.setSelected(params.excludeRound);
            cArea.setSelectedIndex(areaMode);
            fMargin.setValue(margin);
            fBrightness.setValue(display.brightness * 100);
            fContrast.setValue(display.contrast);
            fClarity.setValue(display.clarity);
            fGamma.setValue(display.gamma);
            fClahe.setValue(display.clahe);
            fMedian.setValue(display.median);
            fSmooth.setValue(display.denoise);
            fBackground.setValue(display.background);
            bFlatten.setSelected(display.flatten);
            bPreview.setSelected(preview);
            cStyle.setSelectedIndex(style);
            bLiveMask.setSelected(liveMask);
            bTable.setSelected(outputTable);
            bOutline.setSelected(outputOutline);
            bMask.setSelected(outputMask);
            bAllSlices.setSelected(allSlices);
        } finally {
            updatingControls = false;
        }
    }

    /** Called whenever any control changes. */
    private void controlsChanged() {
        if (updatingControls || dialog == null) return;
        readControls();
        String vk = display.key();
        if (!vk.equals(lastDisplayKey)) {
            lastDisplayKey = vk;
            scheduleDisplay();
        }
        updateAppliedLabel();
        updateAreaLabel();
        refreshDetection();
        rebuildOverlay();
        updatePresetState();
    }

    private String detectionKey() {
        return sigma + "|" + params.key() + "|" + areaMode + "|" + margin + "|" + areaVersion + "|" + editVersion
                + "|" + preview + "|" + (style != STYLE_FILL) + "|" + applied.key();
    }

    private void refreshDetection() {
        String dk = detectionKey();
        if (!dk.equals(lastDetectionKey)) {
            lastDetectionKey = dk;
            scheduleDetection();
        }
    }

    // ---------------------------------------------------------------- buttons

    /** "Apply to detection": detection now runs on the image as currently adjusted. */
    private void applyAdjustments() {
        applied = display;
        updateAppliedLabel();
        refreshDetection();
        updatePresetState();
    }

    /** "Default": adjustment sliders back to neutral, and detection back to the original image. */
    private void resetAdjustments() {
        display = NEUTRAL;
        applied = NEUTRAL;
        syncControls();
        controlsChanged();
    }

    private void addEdit(boolean add) {
        if (imp == null) {
            setStatus("Open an image first");
            return;
        }
        Roi roi = imp.getRoi();
        if (roi == null || !roi.isArea()) {
            setStatus("Draw a selection on the image first (e.g. freehand or oval tool)");
            return;
        }
        edits.add(new Edit(imp.getCurrentSlice(), (Roi) roi.clone(), add));
        editVersion++;
        imp.deleteRoi();
        updateEditsLabel();
        refreshDetection();
    }

    private void undoEdit() {
        if (edits.isEmpty()) return;
        edits.remove(edits.size() - 1);
        editVersion++;
        updateEditsLabel();
        refreshDetection();
    }

    private void clearEdits() {
        if (edits.isEmpty()) return;
        edits.clear();
        editVersion++;
        updateEditsLabel();
        refreshDetection();
    }

    private void setAreaFromSelection() {
        if (imp == null) {
            setStatus("Open an image first");
            return;
        }
        Roi roi = imp.getRoi();
        if (roi == null || !roi.isArea()) {
            setStatus("Draw the analysis area on the image first (e.g. oval tool)");
            return;
        }
        areaRoi = (Roi) roi.clone();
        areaVersion++;
        imp.deleteRoi();
        areaMode = AREA_SELECTION;
        syncControls();
        controlsChanged();
    }

    private void updateAreaLabel() {
        if (areaLabel == null) return;
        boolean showButton = areaMode == AREA_SELECTION;
        if (setAreaButton != null && setAreaButton.isVisible() != showButton) {
            setAreaButton.setVisible(showButton);
            setAreaButton.getParent().revalidate();
        }
        String text;
        if (areaMode == AREA_SELECTION) {
            text = areaRoi == null ? "No area set yet - draw it on the image and press the button above"
                    : "Area set (" + areaRoi.getBounds().width + " x " + areaRoi.getBounds().height + " px, cyan outline)";
        } else if (areaMode == AREA_AUTO) {
            text = "Detected field of view is outlined in cyan";
        } else {
            text = "Analysing the whole image";
        }
        areaLabel.setText(text);
    }

    private void updateEditsLabel() {
        if (editsLabel == null) return;
        int added = 0, removed = 0;
        for (Edit e : edits) if (e.add) added++; else removed++;
        String text = edits.isEmpty() ? "Manual edits: none"
                : "Manual edits: " + added + " added, " + removed + " removed";
        EventQueue.invokeLater(() -> editsLabel.setText(text));
    }

    private void updateAppliedLabel() {
        if (appliedLabel == null) return;
        String text = applied.neutral() ? "Detection uses: original image" : "Detection uses: adjusted image (applied)";
        if (!applied.key().equals(display.key())) text += " - adjustments changed, press Apply to use them";
        final String t = text;
        EventQueue.invokeLater(() -> appliedLabel.setText(t));
    }

    // ---------------------------------------------------------------- presets

    private void loadPreset(String name) {
        Map<String, String> values;
        if (DEFAULT_PRESET.equals(name)) {
            values = defaultSettings();
        } else {
            String stored = loadPresets().get(name);
            if (stored == null) return;
            values = parse(stored);
        }
        fieldsFromMap(values);
        syncControls();
        loadedPresetName = name;
        loadedPresetValues = serialize(currentSettings());
        controlsChanged();
        updateAppliedLabel();
        refreshDetection();
        updatePresetState();
    }

    private void savePreset() {
        String current = serialize(currentSettings());
        if (DEFAULT_PRESET.equals(loadedPresetName)) {
            saveAsNewPreset(current);
            return;
        }
        Object[] choices = {"Update preset", "Save new preset", "Cancel"};
        int r = JOptionPane.showOptionDialog(dialog,
                "Update the preset \"" + loadedPresetName + "\" with the current settings,\nor save them as a new preset?",
                "Save preset", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, choices, choices[0]);
        if (r == 0) {
            Map<String, String> all = loadPresets();
            all.put(loadedPresetName, current);
            if (storePresets(all)) {
                loadedPresetValues = current;
                updatePresetState();
            }
        } else if (r == 1) {
            saveAsNewPreset(current);
        }
    }

    private void saveAsNewPreset(String current) {
        String name = JOptionPane.showInputDialog(dialog, "Name for the new preset:", "Save preset",
                JOptionPane.PLAIN_MESSAGE);
        if (name == null) return;
        name = name.trim();
        if (name.isEmpty()) return;
        if (DEFAULT_PRESET.equalsIgnoreCase(name)) {
            JOptionPane.showMessageDialog(dialog, "\"Default\" is reserved. Please choose another name.");
            return;
        }
        Map<String, String> all = loadPresets();
        if (all.containsKey(name)) {
            int r = JOptionPane.showConfirmDialog(dialog, "A preset called \"" + name + "\" already exists. Replace it?",
                    "Save preset", JOptionPane.YES_NO_OPTION);
            if (r != JOptionPane.YES_OPTION) return;
        }
        all.put(name, current);
        if (!storePresets(all)) return;
        loadedPresetName = name;
        loadedPresetValues = current;
        refreshPresetList(name);
        updatePresetState();
    }

    private void deletePreset() {
        if (DEFAULT_PRESET.equals(loadedPresetName)) {
            JOptionPane.showMessageDialog(dialog, "The Default preset cannot be deleted.");
            return;
        }
        int r = JOptionPane.showConfirmDialog(dialog, "Delete the preset \"" + loadedPresetName + "\"?",
                "Delete preset", JOptionPane.YES_NO_OPTION);
        if (r != JOptionPane.YES_OPTION) return;
        Map<String, String> all = loadPresets();
        all.remove(loadedPresetName);
        if (!storePresets(all)) return;
        loadedPresetName = DEFAULT_PRESET;
        loadedPresetValues = serialize(defaultSettings());
        refreshPresetList(DEFAULT_PRESET);
        updatePresetState();
    }

    private void refreshPresetList(String select) {
        updatingControls = true;
        try {
            cPreset.removeAllItems();
            cPreset.addItem(DEFAULT_PRESET);
            for (String name : new TreeSet<>(loadPresets().keySet())) cPreset.addItem(name);
            cPreset.setSelectedItem(select);
        } finally {
            updatingControls = false;
        }
        refreshBatchPresets();
    }

    private void updatePresetState() {
        if (presetState == null || loadedPresetValues == null) return;
        boolean modified = !serialize(currentSettings()).equals(loadedPresetValues);
        presetState.setText(modified ? "(modified - press Save preset to keep)" : " ");
    }

    /** All settings that a preset stores (not the analysis area outline or manual edits). */
    Map<String, String> currentSettings() {
        Map<String, String> m = new LinkedHashMap<>();
        m.put("sigma", Double.toString(sigma));
        m.put("epsilon", Double.toString(params.epsilon));
        m.put("haloDepth", Integer.toString(params.haloDepth));
        m.put("minSize", Integer.toString(params.minSize));
        m.put("maxHole", Integer.toString(params.maxHole));
        m.put("fillAll", Boolean.toString(params.fillAll));
        m.put("grow", Integer.toString(params.grow));
        m.put("excludeRound", Boolean.toString(params.excludeRound));
        m.put("area", Integer.toString(areaMode));
        m.put("margin", Integer.toString(margin));
        putDisplay(m, "adjust.", display);
        putDisplay(m, "applied.", applied);
        m.put("preview", Boolean.toString(preview));
        m.put("style", Integer.toString(style));
        m.put("liveMask", Boolean.toString(liveMask));
        m.put("table", Boolean.toString(outputTable));
        m.put("outline", Boolean.toString(outputOutline));
        m.put("mask", Boolean.toString(outputMask));
        m.put("allSlices", Boolean.toString(allSlices));
        return m;
    }

    private static void putDisplay(Map<String, String> m, String prefix, Display d) {
        m.put(prefix + "brightness", Double.toString(d.brightness));
        m.put(prefix + "contrast", Double.toString(d.contrast));
        m.put(prefix + "clarity", Double.toString(d.clarity));
        m.put(prefix + "gamma", Double.toString(d.gamma));
        m.put(prefix + "clahe", Double.toString(d.clahe));
        m.put(prefix + "median", Double.toString(d.median));
        m.put(prefix + "smooth", Double.toString(d.denoise));
        m.put(prefix + "background", Double.toString(d.background));
        m.put(prefix + "flatten", Boolean.toString(d.flatten));
    }

    static Map<String, String> defaultSettings() {
        return new PHANTAST_Live().currentSettings();
    }

    /** Sets the settings fields from a stored preset; missing values fall back to the defaults. */
    void fieldsFromMap(Map<String, String> m) {
        Map<String, String> v = defaultSettings();
        v.putAll(m);
        sigma = clamp(num(v, "sigma"), SIGMA_MIN, SIGMA_MAX);
        params = new Params(num(v, "epsilon"), (int) num(v, "haloDepth"), (int) num(v, "minSize"),
                (int) num(v, "maxHole"), bool(v, "fillAll"), (int) num(v, "grow"), bool(v, "excludeRound"));
        areaMode = (int) clamp(num(v, "area"), 0, AREA_NAMES.length - 1);
        margin = (int) clamp(num(v, "margin"), 0, MARGIN_MAX);
        display = getDisplay(v, "adjust.");
        applied = getDisplay(v, "applied.");
        preview = bool(v, "preview");
        style = (int) clamp(num(v, "style"), 0, STYLE_NAMES.length - 1);
        liveMask = bool(v, "liveMask");
        outputTable = bool(v, "table");
        outputOutline = bool(v, "outline");
        outputMask = bool(v, "mask");
        allSlices = bool(v, "allSlices");
    }

    private static Display getDisplay(Map<String, String> v, String prefix) {
        return new Display(num(v, prefix + "brightness"), num(v, prefix + "contrast"), num(v, prefix + "clarity"),
                num(v, prefix + "gamma"), num(v, prefix + "clahe"), num(v, prefix + "median"),
                num(v, prefix + "smooth"), num(v, prefix + "background"), bool(v, prefix + "flatten"));
    }

    private static double num(Map<String, String> v, String key) {
        try {
            return Double.parseDouble(v.get(key));
        } catch (RuntimeException e) {
            return Double.parseDouble(defaultSettings().get(key));
        }
    }

    private static boolean bool(Map<String, String> v, String key) {
        return Boolean.parseBoolean(v.get(key));
    }

    static String serialize(Map<String, String> m) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, String> e : m.entrySet()) sb.append(e.getKey()).append('=').append(e.getValue()).append(';');
        return sb.toString();
    }

    static Map<String, String> parse(String s) {
        Map<String, String> m = new LinkedHashMap<>();
        for (String part : s.split(";")) {
            int eq = part.indexOf('=');
            if (eq > 0) m.put(part.substring(0, eq).trim(), part.substring(eq + 1).trim());
        }
        return m;
    }

    static File presetFile() {
        return new File(Prefs.getPrefsDir(), PRESET_FILE);
    }

    /** Preset name to serialized settings. */
    static Map<String, String> loadPresets() {
        Map<String, String> out = new LinkedHashMap<>();
        File f = presetFile();
        if (!f.exists()) return out;
        Properties p = new Properties();
        try (InputStream in = new FileInputStream(f)) {
            p.load(in);
        } catch (Exception e) {
            IJ.log("PHANTAST Live: could not read presets from " + f + ": " + e);
            return out;
        }
        for (String name : p.stringPropertyNames()) out.put(name, p.getProperty(name));
        return out;
    }

    static boolean storePresets(Map<String, String> presets) {
        Properties p = new Properties();
        p.putAll(presets);
        File f = presetFile();
        try {
            File dir = f.getParentFile();
            if (dir != null && !dir.exists()) dir.mkdirs();
            try (OutputStream out = new FileOutputStream(f)) {
                p.store(out, "PHANTAST Live presets");
            }
            return true;
        } catch (Exception e) {
            IJ.error("PHANTAST Live", "Could not save presets to " + f + ":\n" + e.getMessage());
            return false;
        }
    }

    // ---------------------------------------------------------------- live preview

    private void scheduleDisplay() {
        final int g = displayGen.incrementAndGet();
        final ImagePlus target = imp;
        if (target == null) {
            layerView = null;
            return;
        }
        final Display d = display;
        final int slice = target.getCurrentSlice();
        final BooleanSupplier stale = () -> displayGen.get() != g;
        displayWorker.submit(() -> {
            if (stale.getAsBoolean()) return;
            try {
                ImageRoi view = renderDisplay(target, slice, d, stale);
                if (stale.getAsBoolean()) return;
                layerView = view;
                EventQueue.invokeLater(this::rebuildOverlay);
            } catch (Cancelled x) {
                // superseded
            } catch (OutOfMemoryError oom) {
                vRaw = null; vBase = null; vBaseBlur = null; vBaseKey = null; vSlice = -1;
                setStatus("Not enough memory for image adjustments");
            } catch (Throwable t) {
                IJ.log("PHANTAST Live display error: " + t);
            }
        });
    }

    private void scheduleDetection() {
        final int g = detectGen.incrementAndGet();
        final ImagePlus target = imp;
        if (target == null) {
            layerFill = null; layerOutline = null; layerArea = null; layerMask = null;
            setStatus("Confluency: open an image to see it");
            return;
        }
        if (!preview) {
            layerFill = null; layerOutline = null; layerArea = null; layerMask = null;
            setStatus("Confluency: (turn on Live preview)");
            EventQueue.invokeLater(this::rebuildOverlay);
            return;
        }
        setStatus("Confluency: calculating...");
        final double s = sigma;
        final Params p = params;
        final boolean wantOutline = style != STYLE_FILL;
        final Display adj = applied;
        final int mode = areaMode, mg = margin;
        final Roi area = areaRoi;
        final int slice = target.getCurrentSlice();
        final List<Edit> sliceEdits = editsFor(slice);
        final boolean stack = target.getStackSize() > 1;
        final BooleanSupplier stale = () -> detectGen.get() != g;
        detectWorker.submit(() -> {
            try {
                Thread.sleep(40); // let slider drags settle
            } catch (InterruptedException ie) {
                return;
            }
            if (stale.getAsBoolean()) return;
            try {
                Detection det = detectPreview(target, slice, adj, s, p, mode, mg, area, sliceEdits, stale);
                if (det.area.error != null) {
                    if (stale.getAsBoolean()) return;
                    layerFill = null; layerOutline = null; layerArea = null; layerMask = null;
                    setStatus(det.area.error);
                    EventQueue.invokeLater(this::rebuildOverlay);
                    return;
                }
                ImageRoi fill = tint(det.mask);
                Roi outline = null;
                if (wantOutline) {
                    outline = selection(det.mask);
                    if (outline != null) outline.setStrokeColor(CELL_OUTLINE);
                }
                Roi areaOutline = null;
                if (det.area.outline != null) {
                    areaOutline = (Roi) det.area.outline.clone();
                    areaOutline.setStrokeColor(AREA_OUTLINE);
                    areaOutline.setStrokeWidth(2);
                }
                if (stack) {
                    fill.setPosition(slice);
                    if (outline != null) outline.setPosition(slice);
                    if (areaOutline != null) areaOutline.setPosition(slice);
                }
                if (stale.getAsBoolean()) return;
                layerFill = fill;
                layerOutline = outline;
                layerArea = areaOutline;
                layerMask = det.mask;
                final double c = det.confluency;
                EventQueue.invokeLater(() -> {
                    if (stale.getAsBoolean()) return;
                    status.setText(String.format("Confluency: %.1f %%", c));
                    rebuildOverlay();
                });
            } catch (Cancelled x) {
                // superseded by newer settings
            } catch (OutOfMemoryError oom) {
                dImage = null; dInput = null; dContrast = null; dDirections = null; dAutoArea = null; dSlice = -1;
                setStatus("Not enough memory for preview");
            } catch (Throwable t) {
                setStatus("Preview failed: " + t.getClass().getSimpleName());
                IJ.log("PHANTAST Live preview error: " + t);
            }
        });
    }

    private List<Edit> editsFor(int slice) {
        List<Edit> out = new ArrayList<>();
        boolean single = imp == null || imp.getStackSize() == 1;
        for (Edit e : edits) if (e.slice == slice || single) out.add(e);
        return out;
    }

    /** Detection for the preview, reusing cached steps when only some settings changed. */
    private Detection detectPreview(ImagePlus target, int slice, Display adj, double s, Params p, int mode, int mg,
                                    Roi areaSel, List<Edit> sliceEdits, BooleanSupplier stale) {
        if (target != dImp || slice != dSlice || dImage == null) {
            dImage = toFloat(target.getStack().getProcessor(slice));
            dImp = target;
            dSlice = slice;
            dAdjustKey = null;
            dInput = null;
            dAutoArea = null;
        }
        if (dInput == null || !adj.key().equals(dAdjustKey)) {
            dInput = detectionInput(target, target.getStack().getProcessor(slice), adj);
            dAdjustKey = adj.key();
            dSigma = Double.NaN;
            dContrast = null;
            dDirections = null;
        }
        int w = dImage.getWidth(), h = dImage.getHeight();
        Area area;
        if (mode == AREA_AUTO) {
            if (dAutoArea == null || dAutoMargin != mg) {
                dAutoArea = autoArea(dImage, mg);
                dAutoMargin = mg;
            }
            area = dAutoArea;
        } else {
            area = buildArea(mode, areaSel, dImage, mg);
        }
        if (area.error != null) return new Detection(null, 0, area);
        if (stale.getAsBoolean()) throw new Cancelled();
        if (dContrast == null || dSigma != s) {
            dContrast = localContrast(dInput, s);
            dSigma = s;
        }
        if (stale.getAsBoolean()) throw new Cancelled();
        if (p.haloDepth > 0 && dDirections == null) dDirections = directions(dInput);
        float[] raw = (float[]) dImage.getPixels();
        ByteProcessor mask = segment(dContrast, dDirections, raw, w, h, p, area.mask, sliceEdits, stale);
        return new Detection(mask, confluency(mask, area.mask), area);
    }

    /** Draws all current layers on the image; runs on the event thread. */
    private void rebuildOverlay() {
        if (imp == null) {
            closeLiveMask();
            return;
        }
        Overlay ov = originalOverlay == null ? new Overlay() : originalOverlay.duplicate();
        ImageRoi view = layerView;
        if (view != null) {
            if (imp.getStackSize() > 1) view.setPosition(imp.getCurrentSlice());
            ov.add(view);
        }
        if (preview && !(hideDetection && onAdjustTab)) {
            ImageRoi fill = layerFill;
            Roi outline = layerOutline, area = layerArea;
            if (fill != null && style != STYLE_OUTLINE) ov.add(fill);
            if (outline != null && style != STYLE_FILL) ov.add(outline);
            if (area != null) ov.add(area);
        }
        imp.setOverlay(ov.size() == 0 ? originalOverlay : ov);
        updateLiveMask();
    }

    private void updateLiveMask() {
        ByteProcessor mask = layerMask;
        if (!liveMask || !preview || mask == null) {
            closeLiveMask();
            return;
        }
        if (liveMaskImp == null || liveMaskImp.getWindow() == null) {
            liveMaskImp = new ImagePlus(imp.getShortTitle() + " - live mask", mask.duplicate());
            liveMaskImp.show();
            if (dialog != null) dialog.toFront();
        } else {
            liveMaskImp.setProcessor(mask.duplicate());
            liveMaskImp.updateAndDraw();
        }
    }

    private void closeLiveMask() {
        if (liveMaskImp == null) return;
        ImagePlus m = liveMaskImp;
        liveMaskImp = null;
        m.changes = false;
        if (m.getWindow() != null) m.close();
    }

    /** Builds the adjusted greyscale view; only the per-pixel pass runs when a fast slider moves. */
    private ImageRoi renderDisplay(ImagePlus target, int slice, Display d, BooleanSupplier stale) {
        if (d.neutral()) return null;
        if (target != vImp || slice != vSlice || vRaw == null) {
            vRaw = displayFloat(target, target.getStack().getProcessor(slice));
            vImp = target;
            vSlice = slice;
            vBaseKey = null;
        }
        if (!d.baseKey().equals(vBaseKey)) {
            FloatProcessor base = d.buildBase(vRaw);
            if (stale.getAsBoolean()) throw new Cancelled();
            vBase = (float[]) base.getPixels();
            vBaseBlur = null;
            vMean = mean(vBase);
            vBaseKey = d.baseKey();
        }
        if (d.clarity > 0 && vBaseBlur == null) {
            FloatProcessor blur = new FloatProcessor(vRaw.getWidth(), vRaw.getHeight(), vBase.clone());
            new GaussianBlur().blurGaussian(blur, Display.CLARITY_RADIUS, Display.CLARITY_RADIUS, 0.01);
            vBaseBlur = (float[]) blur.getPixels();
        }
        if (stale.getAsBoolean()) throw new Cancelled();
        byte[] out = renderPixels(vBase, d.clarity > 0 ? vBaseBlur : null, vMean, d);
        ImageRoi roi = new ImageRoi(0, 0, new ByteProcessor(vRaw.getWidth(), vRaw.getHeight(), out));
        roi.setOpacity(1.0);
        return roi;
    }

    private void setStatus(String text) {
        if (status == null) return;
        EventQueue.invokeLater(() -> status.setText(text));
    }

    private static ImageRoi tint(ByteProcessor mask) {
        byte[] m = (byte[]) mask.getPixels();
        int[] rgb = new int[m.length];
        for (int i = 0; i < m.length; i++) if (m[i] != 0) rgb[i] = 0x00E676;
        ImageRoi roi = new ImageRoi(0, 0, new ColorProcessor(mask.getWidth(), mask.getHeight(), rgb));
        roi.setZeroTransparent(true);
        roi.setOpacity(0.4);
        return roi;
    }

    // ---------------------------------------------------------------- image adjustments

    /** Brightness, contrast (around the mean), clarity (unsharp mask) and gamma in one pass. */
    static byte[] renderPixels(float[] base, float[] baseBlur, double mean, Display d) {
        byte[] lut = new byte[256];
        for (int i = 0; i < 256; i++) lut[i] = (byte) Math.round(255 * Math.pow(i / 255.0, d.gamma));
        byte[] out = new byte[base.length];
        double c = d.contrast;
        double offset = mean * (1 - c) + d.brightness;
        double k = d.clarity * c;
        for (int i = 0; i < base.length; i++) {
            double v = base[i] * c + offset;
            if (baseBlur != null) v += k * (base[i] - baseBlur[i]);
            int g = (int) Math.round(v * 255);
            out[i] = lut[g < 0 ? 0 : g > 255 ? 255 : g];
        }
        return out;
    }

    /** The image as currently displayed (respecting its display range), scaled to 0..1. */
    static FloatProcessor displayFloat(ImagePlus imp, ImageProcessor ip) {
        ByteProcessor bp;
        if (ip instanceof ByteProcessor) {
            bp = (ByteProcessor) ip;
        } else if (ip instanceof ColorProcessor) {
            bp = (ByteProcessor) ip.convertToByte(false);
        } else {
            ImageProcessor copy = ip.duplicate();
            copy.setMinAndMax(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());
            bp = (ByteProcessor) copy.convertToByte(true);
        }
        byte[] b = (byte[]) bp.getPixels();
        float[] f = new float[b.length];
        for (int i = 0; i < b.length; i++) f[i] = (b[i] & 0xff) / 255f;
        return new FloatProcessor(bp.getWidth(), bp.getHeight(), f);
    }

    /** Input image for detection: original pixels, or the adjusted view when adjustments are applied. */
    static FloatProcessor detectionInput(ImagePlus imp, ImageProcessor ip, Display adj) {
        if (adj.neutral()) return toFloat(ip);
        return adj.apply(displayFloat(imp, ip));
    }

    /** Divides out a heavily blurred background to correct uneven illumination. */
    static FloatProcessor flatten(FloatProcessor img) {
        int w = img.getWidth(), h = img.getHeight();
        FloatProcessor bg = (FloatProcessor) img.duplicate();
        double s = Math.max(w, h) / 20.0;
        new GaussianBlur().blurGaussian(bg, s, s, 0.01);
        return divideByBackground(img, (float[]) bg.getPixels());
    }

    /** Rolling-ball background (light background) divided out, keeping the image brightness. */
    static FloatProcessor rollingBallFlatten(FloatProcessor img, double radius) {
        FloatProcessor bg = (FloatProcessor) img.duplicate();
        new BackgroundSubtracter().rollingBallBackground(bg, radius, true, true, false, true, true);
        return divideByBackground(img, (float[]) bg.getPixels());
    }

    static FloatProcessor divideByBackground(FloatProcessor img, float[] b) {
        float[] v = ((float[]) img.getPixels()).clone();
        double mean = mean(v);
        double floor = Math.max(1e-6, 0.05 * mean);
        for (int i = 0; i < v.length; i++) v[i] = (float) (v[i] / Math.max(b[i], floor) * mean);
        return new FloatProcessor(img.getWidth(), img.getHeight(), v);
    }

    /**
     * Contrast-limited adaptive histogram equalisation. {@code slope} is the clip limit
     * (higher = stronger local contrast); tiles are blended bilinearly.
     */
    static FloatProcessor clahe(FloatProcessor img, double slope, int tile) {
        int w = img.getWidth(), h = img.getHeight();
        float[] v = (float[]) img.getPixels();
        int bins = 256;
        int tx = Math.max(1, (w + tile - 1) / tile), ty = Math.max(1, (h + tile - 1) / tile);
        float[][] maps = new float[tx * ty][];
        int[] bin = new int[v.length];
        for (int i = 0; i < v.length; i++) {
            int b = (int) (v[i] * (bins - 1) + 0.5f);
            bin[i] = b < 0 ? 0 : b >= bins ? bins - 1 : b;
        }
        for (int j = 0; j < ty; j++) {
            for (int i = 0; i < tx; i++) {
                int x0 = i * w / tx, x1 = (i + 1) * w / tx, y0 = j * h / ty, y1 = (j + 1) * h / ty;
                double[] hist = new double[bins];
                int count = 0;
                for (int y = y0; y < y1; y++)
                    for (int x = x0; x < x1; x++) { hist[bin[y * w + x]]++; count++; }
                double limit = Math.max(1, 1 + slope * count / (double) bins);
                double excess = 0;
                for (int b = 0; b < bins; b++) {
                    if (hist[b] > limit) { excess += hist[b] - limit; hist[b] = limit; }
                }
                double add = excess / bins;
                float[] map = new float[bins];
                double cdf = 0;
                for (int b = 0; b < bins; b++) {
                    cdf += hist[b] + add;
                    map[b] = (float) (cdf / Math.max(1, count));
                }
                maps[j * tx + i] = map;
            }
        }
        float[] out = new float[v.length];
        double cw = w / (double) tx, chh = h / (double) ty;
        for (int y = 0; y < h; y++) {
            double fy = (y + 0.5) / chh - 0.5;
            int j0 = (int) Math.floor(fy);
            double wy = fy - j0;
            int ja = Math.max(0, Math.min(ty - 1, j0)), jb = Math.max(0, Math.min(ty - 1, j0 + 1));
            for (int x = 0; x < w; x++) {
                double fx = (x + 0.5) / cw - 0.5;
                int i0 = (int) Math.floor(fx);
                double wx = fx - i0;
                int ia = Math.max(0, Math.min(tx - 1, i0)), ib = Math.max(0, Math.min(tx - 1, i0 + 1));
                int b = bin[y * w + x];
                double top = maps[ja * tx + ia][b] * (1 - wx) + maps[ja * tx + ib][b] * wx;
                double bottom = maps[jb * tx + ia][b] * (1 - wx) + maps[jb * tx + ib][b] * wx;
                out[y * w + x] = (float) (top * (1 - wy) + bottom * wy);
            }
        }
        return new FloatProcessor(w, h, out);
    }

    static double mean(float[] v) {
        double s = 0;
        for (float f : v) s += f;
        return v.length == 0 ? 0 : s / v.length;
    }

    // ---------------------------------------------------------------- final run

    private void finalRun() {
        int n = imp.getStackSize();
        int[] slices;
        if (n > 1 && allSlices) {
            slices = new int[n];
            for (int i = 0; i < n; i++) slices[i] = i + 1;
        } else {
            slices = new int[] {imp.getCurrentSlice()};
        }
        if (areaMode == AREA_SELECTION && (areaRoi == null || !areaRoi.isArea())) {
            IJ.error("PHANTAST Live", "Analysis area is set to 'Selection', but no area has been set.\n"
                    + "Draw the area on the image and press 'Set area from selection'.");
            return;
        }
        ResultsTable rt = outputTable ? ResultsTable.getResultsTable() : null;
        Overlay ov = null;
        if (outputOutline) ov = originalOverlay == null ? new Overlay() : originalOverlay.duplicate();
        ImageStack maskStack = outputMask ? new ImageStack(imp.getWidth(), imp.getHeight()) : null;
        double last = 0;

        for (int i = 0; i < slices.length; i++) {
            int s = slices[i];
            IJ.showStatus("PHANTAST: image " + (i + 1) + "/" + slices.length);
            IJ.showProgress(i, slices.length);
            List<Edit> sliceEdits = editsFor(s);
            Measurement r = measure(imp, s, areaRoi, sliceEdits);
            if (r.area.error != null) {
                IJ.log("PHANTAST Live, slice " + s + ": " + r.area.error);
                continue;
            }
            ByteProcessor mask = r.mask;
            Area area = r.area;
            last = r.confluency;
            if (rt != null) {
                rt.incrementCounter();
                rt.addValue("Image", imp.getTitle());
                if (n > 1) rt.addValue("Slice", s);
                addSettingsColumns(rt, r, sliceEdits.size());
            }
            if (ov != null) {
                Roi roi = selection(mask);
                if (roi != null) {
                    roi.setStrokeColor(CELL_OUTLINE);
                    if (n > 1) roi.setPosition(s);
                    ov.add(roi);
                }
                if (area.outline != null) {
                    Roi edge = (Roi) area.outline.clone();
                    edge.setStrokeColor(AREA_OUTLINE);
                    edge.setStrokeWidth(2);
                    if (n > 1) edge.setPosition(s);
                    ov.add(edge);
                }
            }
            if (maskStack != null) maskStack.addSlice(imp.getStack().getShortSliceLabel(s), mask);
        }
        IJ.showProgress(1.0);

        if (rt != null) rt.show("Results");
        if (ov != null) imp.setOverlay(ov);
        if (maskStack != null && maskStack.getSize() > 0)
            new ImagePlus(imp.getShortTitle() + " - mask", maskStack).show();
        IJ.showStatus(String.format("PHANTAST: confluency %.1f %%", last));
    }

    // ---------------------------------------------------------------- measurement (shared by OK and batch)

    /** Detection result for one image or slice. */
    static final class Measurement {
        final ByteProcessor mask;
        final double confluency;
        final Area area;
        final long areaPixels;
        Measurement(ByteProcessor mask, double confluency, Area area, long areaPixels) {
            this.mask = mask;
            this.confluency = confluency;
            this.area = area;
            this.areaPixels = areaPixels;
        }
    }

    /** Runs the full detection with this instance's settings on one slice of an image. */
    Measurement measure(ImagePlus im, int slice, Roi areaSel, List<Edit> sliceEdits) {
        ImageProcessor ip = im.getStack().getProcessor(slice);
        FloatProcessor img = toFloat(ip);
        int w = img.getWidth(), h = img.getHeight();
        Area area = buildArea(areaMode, areaSel, img, margin);
        if (area.error != null) return new Measurement(null, 0, area, 0);
        FloatProcessor input = detectionInput(im, ip, applied);
        float[] lc = localContrast(input, sigma);
        byte[] dir = params.haloDepth > 0 ? directions(input) : null;
        ByteProcessor mask = segment(lc, dir, (float[]) img.getPixels(), w, h, params, area.mask, sliceEdits, () -> false);
        long areaPixels = area.mask == null ? (long) w * h : count(area.mask);
        return new Measurement(mask, confluency(mask, area.mask), area, areaPixels);
    }

    /** Confluency and every setting used, so each row can be reported and reproduced. */
    void addSettingsColumns(ResultsTable rt, Measurement r, int manualEdits) {
        rt.addValue("Confluency (%)", r.confluency);
        rt.addValue("Analysis area", areaMode == AREA_WHOLE ? "Whole image"
                : areaMode == AREA_SELECTION ? "Selection" : "Auto circle");
        rt.addValue("Area (px)", r.areaPixels);
        rt.addValue("Sigma", sigma);
        rt.addValue("Epsilon", params.epsilon);
        rt.addValue("Halo correction", params.haloDepth == 0 ? "Off"
                : params.haloDepth >= HALO_FULL ? "Full" : params.haloDepth + " px");
        rt.addValue("Min cell size (px)", params.minSize);
        rt.addValue("Fill holes", params.fillAll ? "All" : params.maxHole + " px");
        rt.addValue("Grow/shrink (px)", params.grow);
        rt.addValue("Exclude round bright", params.excludeRound ? "Yes" : "No");
        rt.addValue("Manual edits", manualEdits);
        rt.addValue("Image adjustments", applied.describe());
    }

    /** The image as displayed, with the cell outline in yellow and the analysis area in cyan. */
    static ColorProcessor outlineImage(ImagePlus im, int slice, Measurement r) {
        FloatProcessor view = displayFloat(im, im.getStack().getProcessor(slice));
        float[] v = (float[]) view.getPixels();
        byte[] b = new byte[v.length];
        for (int i = 0; i < v.length; i++) b[i] = (byte) Math.round(v[i] * 255);
        ColorProcessor cp = (ColorProcessor) new ByteProcessor(view.getWidth(), view.getHeight(), b).convertToRGB();
        Roi cells = selection(r.mask);
        int width = Math.max(1, Math.max(cp.getWidth(), cp.getHeight()) / 1000);
        if (cells != null) {
            cp.setColor(CELL_OUTLINE);
            cp.setLineWidth(width);
            cp.draw(cells);
        }
        if (r.area.outline != null) {
            cp.setColor(AREA_OUTLINE);
            cp.setLineWidth(width + 1);
            cp.draw(r.area.outline);
        }
        return cp;
    }

    // ---------------------------------------------------------------- batch processing

    static final String[] IMAGE_EXTENSIONS = {".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp", ".gif"};
    static final String RESULTS_PREFIX = "PHANTAST_Live_results_";

    interface BatchProgress {
        void update(int done, int total, String message);
    }

    static final class BatchResult {
        File outputDir;
        ResultsTable table;
        int images, failed;
        boolean cancelled;
        String error;
    }

    /** Image files in a folder (optionally its subfolders), skipping earlier result folders. */
    static List<File> findImages(File dir, boolean recurse) {
        List<File> out = new ArrayList<>();
        File[] entries = dir.listFiles();
        if (entries == null) return out;
        Arrays.sort(entries);
        for (File f : entries) {
            if (f.isDirectory()) {
                if (recurse && !f.getName().startsWith(RESULTS_PREFIX)) out.addAll(findImages(f, true));
            } else {
                String name = f.getName().toLowerCase();
                if (name.startsWith(".")) continue;
                for (String ext : IMAGE_EXTENSIONS) {
                    if (name.endsWith(ext)) {
                        out.add(f);
                        break;
                    }
                }
            }
        }
        return out;
    }

    /**
     * Analyses every image in {@code dir} with the given settings and saves a results CSV
     * (plus optional masks and outline images) in a new results folder inside {@code dir}.
     */
    static BatchResult runBatch(File dir, boolean recurse, Map<String, String> settings, Roi areaSel,
                                boolean saveMasks, boolean saveOutlines, BooleanSupplier cancelled,
                                BatchProgress progress) {
        BatchResult res = new BatchResult();
        PHANTAST_Live b = new PHANTAST_Live();
        b.fieldsFromMap(settings);
        if (b.areaMode == AREA_SELECTION && (areaSel == null || !areaSel.isArea())) {
            res.error = "These settings use 'Selection drawn on image' but no area is set. "
                    + "Set the area in tab 1, or use Whole image / Auto-detect.";
            return res;
        }
        List<File> files = findImages(dir, recurse);
        if (files.isEmpty()) {
            res.error = "No images (" + String.join(", ", IMAGE_EXTENSIONS) + ") found in " + dir;
            return res;
        }
        String stamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new java.util.Date());
        File out = new File(dir, RESULTS_PREFIX + stamp);
        File maskDir = new File(out, "masks"), outlineDir = new File(out, "outlines");
        if (!out.mkdirs() || (saveMasks && !maskDir.mkdirs()) || (saveOutlines && !outlineDir.mkdirs())) {
            res.error = "Could not create the results folder " + out;
            return res;
        }
        ResultsTable rt = new ResultsTable();
        rt.showRowNumbers(false);
        rt.setNaNEmptyCells(true); // rows for unreadable files must not look like 0 % confluency
        String base = dir.getAbsolutePath();

        for (int i = 0; i < files.size(); i++) {
            if (cancelled.getAsBoolean()) {
                res.cancelled = true;
                break;
            }
            File f = files.get(i);
            String rel = f.getAbsolutePath().substring(base.length()).replaceFirst("^[\\\\/]+", "");
            progress.update(i, files.size(), rel);
            try {
                ImagePlus im = IJ.openImage(f.getPath());
                if (im == null) throw new IllegalStateException("could not open this file as an image");
                int n = im.getStackSize();
                String stem = rel.replaceAll("[\\\\/]", "_").replaceFirst("\\.[^.]+$", "");
                for (int s = 1; s <= n; s++) {
                    Measurement r = b.measure(im, s, areaSel, new ArrayList<>());
                    rt.incrementCounter();
                    rt.addValue("File", rel);
                    rt.addValue("Slice", s);
                    if (r.area.error != null) {
                        rt.addValue("Note", r.area.error);
                        continue;
                    }
                    b.addSettingsColumns(rt, r, 0);
                    rt.addValue("Note", "");
                    String name = n > 1 ? stem + "_slice" + s : stem;
                    if (saveMasks)
                        new ij.io.FileSaver(new ImagePlus(name, r.mask)).saveAsPng(new File(maskDir, name + "_mask.png").getPath());
                    if (saveOutlines)
                        new ij.io.FileSaver(new ImagePlus(name, outlineImage(im, s, r)))
                                .saveAsPng(new File(outlineDir, name + "_outline.png").getPath());
                }
                im.flush();
                res.images++;
            } catch (Throwable t) {
                res.failed++;
                rt.incrementCounter();
                rt.addValue("File", rel);
                rt.addValue("Note", "Error: " + t.getMessage());
            }
        }
        progress.update(res.images + res.failed, files.size(), res.cancelled ? "Stopped" : "Saving results");

        try {
            rt.saveAs(new File(out, "confluency_results.csv").getPath());
        } catch (java.io.IOException e) {
            res.error = "Could not save confluency_results.csv: " + e.getMessage();
        }
        List<String> lines = new ArrayList<>();
        lines.add("PHANTAST Live batch settings");
        lines.add("Folder: " + dir.getAbsolutePath());
        lines.add("Images analysed: " + res.images + ", failed: " + res.failed + (res.cancelled ? " (stopped early)" : ""));
        for (Map.Entry<String, String> e : b.currentSettings().entrySet()) lines.add(e.getKey() + " = " + e.getValue());
        try {
            java.nio.file.Files.write(new File(out, "settings_used.txt").toPath(), lines,
                    java.nio.charset.StandardCharsets.UTF_8);
        } catch (Exception e) {
            IJ.log("PHANTAST Live: could not write settings_used.txt: " + e);
        }
        res.outputDir = out;
        res.table = rt;
        return res;
    }

    static Roi selection(ByteProcessor mask) {
        ByteProcessor copy = (ByteProcessor) mask.duplicate();
        copy.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
        return new ThresholdToSelection().convert(copy);
    }

    // ---------------------------------------------------------------- analysis area

    static Area buildArea(int mode, Roi selection, FloatProcessor img, int margin) {
        int w = img.getWidth(), h = img.getHeight();
        if (mode == AREA_SELECTION) {
            if (selection == null || !selection.isArea())
                return new Area(null, null, "Draw the area on the image, then press 'Set area from selection'");
            byte[] m = roiMask(selection, w, h);
            if (count(m) == 0) return new Area(null, null, "The selected area is outside the image");
            return new Area(m, selection, null);
        }
        if (mode == AREA_AUTO) return autoArea(img, margin);
        return new Area(null, null, null);
    }

    static byte[] roiMask(Roi roi, int w, int h) {
        ByteProcessor bp = new ByteProcessor(w, h);
        bp.setColor(255);
        bp.fill(roi);
        return (byte[]) bp.getPixels();
    }

    /**
     * Finds the bright circular field of view in images captured through the eyepiece:
     * everything clearly brighter than the dark surround, largest piece, holes filled,
     * then shrunk by {@code margin} pixels so the bright/dark edge is not counted as cells.
     */
    static Area autoArea(FloatProcessor img, int margin) {
        int w = img.getWidth(), h = img.getHeight(), n = w * h;
        FloatProcessor blurred = (FloatProcessor) img.duplicate();
        double s = Math.max(2, Math.max(w, h) / 200.0);
        new GaussianBlur().blurGaussian(blurred, s, s, 0.01);
        float[] p = (float[]) blurred.getPixels();

        // Reference brightness from the image centre, which is inside the field of view
        double sum = 0;
        int cnt = 0;
        for (int y = h / 3; y < 2 * h / 3; y++)
            for (int x = w / 3; x < 2 * w / 3; x++) { sum += p[y * w + x]; cnt++; }
        double threshold = 0.4 * sum / Math.max(1, cnt);

        byte[] m = new byte[n];
        for (int i = 0; i < n; i++) m[i] = p[i] > threshold ? FG : 0;
        keepLargestRegion(m, w, h);
        flipSmallRegions(m, w, h, (byte) 0, FG, Integer.MAX_VALUE, false, true); // fill all holes

        if (margin > 0) {
            ByteProcessor bp = new ByteProcessor(w, h, m);
            FloatProcessor edm = new EDM().makeFloatEDM(bp, 0, false);
            float[] d = (float[]) edm.getPixels();
            for (int i = 0; i < n; i++) if (d[i] <= margin) m[i] = 0;
        }
        long area = count(m);
        if (area < 0.05 * n) return new Area(null, null, "Could not find a bright field of view - try 'Selection'");
        if (area >= n) return new Area(null, null, null); // no dark surround: whole image

        ByteProcessor bp = new ByteProcessor(w, h, m.clone());
        return new Area(m, selection(bp), null);
    }

    static void keepLargestRegion(byte[] m, int w, int h) {
        int n = w * h;
        int[] label = new int[n];
        int[] stack = new int[n];
        int next = 0, best = 0, bestSize = 0;
        for (int start = 0; start < n; start++) {
            if (label[start] != 0 || m[start] != FG) continue;
            next++;
            int sp = 0, size = 0;
            stack[sp++] = start;
            label[start] = next;
            while (sp > 0) {
                int p = stack[--sp];
                size++;
                int x = p % w, y = p / w;
                if (x > 0 && label[p - 1] == 0 && m[p - 1] == FG) { label[p - 1] = next; stack[sp++] = p - 1; }
                if (x < w - 1 && label[p + 1] == 0 && m[p + 1] == FG) { label[p + 1] = next; stack[sp++] = p + 1; }
                if (y > 0 && label[p - w] == 0 && m[p - w] == FG) { label[p - w] = next; stack[sp++] = p - w; }
                if (y < h - 1 && label[p + w] == 0 && m[p + w] == FG) { label[p + w] = next; stack[sp++] = p + w; }
            }
            if (size > bestSize) { bestSize = size; best = next; }
        }
        for (int i = 0; i < n; i++) if (label[i] != best) m[i] = 0;
    }

    // ---------------------------------------------------------------- detection algorithm

    /** Converts any image type to a float copy scaled to 0..1. */
    static FloatProcessor toFloat(ImageProcessor ip) {
        if (ip instanceof ColorProcessor) ip = ip.convertToByte(false);
        FloatProcessor fp = (FloatProcessor) ip.convertToFloat().duplicate();
        float[] px = (float[]) fp.getPixels();
        float max = 0;
        for (float v : px) if (v > max) max = v;
        if (max > 0) for (int i = 0; i < px.length; i++) px[i] /= max;
        return fp;
    }

    /** Local contrast: local standard deviation divided by local mean (Gaussian window). */
    static float[] localContrast(FloatProcessor img, double sigma) {
        int n = img.getWidth() * img.getHeight();
        float[] px = (float[]) img.getPixels();
        FloatProcessor mean = (FloatProcessor) img.duplicate();
        FloatProcessor sq = new FloatProcessor(img.getWidth(), img.getHeight());
        float[] s = (float[]) sq.getPixels();
        for (int i = 0; i < n; i++) s[i] = px[i] * px[i];
        GaussianBlur gb = new GaussianBlur();
        gb.blurGaussian(mean, sigma, sigma, 0.002);
        gb.blurGaussian(sq, sigma, sigma, 0.002);
        float[] m = (float[]) mean.getPixels();
        float[] lc = new float[n];
        for (int i = 0; i < n; i++) {
            double mu = m[i];
            if (mu > 0) {
                double var = s[i] - mu * mu;
                lc[i] = var > 0 ? (float) (Math.sqrt(var) / mu) : 0f;
            }
        }
        return lc;
    }

    /** Direction of the strongest Kirsch compass response at each pixel (0..7). */
    static byte[] directions(FloatProcessor img) {
        int n = img.getWidth() * img.getHeight();
        float[] best = new float[n];
        byte[] dir = new byte[n];
        for (int k = 0; k < 8; k++) {
            FloatProcessor c = (FloatProcessor) img.duplicate();
            c.convolve(KIRSCH[k], 3, 3);
            float[] r = (float[]) c.getPixels();
            for (int i = 0; i < n; i++) {
                if (k == 0 || r[i] > best[i]) {
                    best[i] = r[i];
                    dir[i] = (byte) k;
                }
            }
        }
        return dir;
    }

    /**
     * Full detection after the local-contrast step: threshold, remove small objects, fill holes,
     * halo correction, grow/shrink, exclude round bright cells, manual edits, then limit to the area.
     */
    static ByteProcessor segment(float[] lc, byte[] dir, float[] raw, int w, int h, Params p, byte[] area,
                                 List<Edit> sliceEdits, BooleanSupplier stale) {
        byte[] m = new byte[w * h];
        for (int i = 0; i < m.length; i++) m[i] = lc[i] > p.epsilon && (area == null || area[i] != 0) ? FG : 0;
        if (stale.getAsBoolean()) throw new Cancelled();
        if (p.minSize > 0) flipSmallRegions(m, w, h, FG, (byte) 0, p.minSize, true, false);
        if (stale.getAsBoolean()) throw new Cancelled();
        int holeLimit = p.fillAll ? Integer.MAX_VALUE : p.maxHole;
        if (holeLimit > 0) flipSmallRegions(m, w, h, (byte) 0, FG, holeLimit, false, true);
        if (p.haloDepth > 0 && dir != null)
            haloCorrection(m, dir, w, h, p.haloDepth >= HALO_FULL ? Integer.MAX_VALUE : p.haloDepth, stale);
        if (stale.getAsBoolean()) throw new Cancelled();
        if (p.grow != 0) growShrink(m, w, h, p.grow);
        if (p.excludeRound && raw != null) excludeRoundBright(m, raw, w, h);
        if (sliceEdits != null) applyEdits(m, w, h, sliceEdits);
        if (area != null) for (int i = 0; i < m.length; i++) if (area[i] == 0) m[i] = 0;
        return new ByteProcessor(w, h, m);
    }

    /** Positive r grows every cell outline by r pixels; negative r shrinks it. */
    static void growShrink(byte[] m, int w, int h, int r) {
        if (r > 0) {
            byte[] inv = new byte[m.length];
            for (int i = 0; i < m.length; i++) inv[i] = m[i] == 0 ? FG : 0;
            float[] d = (float[]) new EDM().makeFloatEDM(new ByteProcessor(w, h, inv), 0, false).getPixels();
            for (int i = 0; i < m.length; i++) if (m[i] == 0 && d[i] <= r) m[i] = FG;
        } else if (r < 0) {
            float[] d = (float[]) new EDM().makeFloatEDM(new ByteProcessor(w, h, m.clone()), 0, false).getPixels();
            for (int i = 0; i < m.length; i++) if (m[i] != 0 && d[i] <= -r) m[i] = 0;
        }
    }

    /** Applies manual corrections in order: add sets the region to cells, remove clears it. */
    static void applyEdits(byte[] m, int w, int h, List<Edit> sliceEdits) {
        for (Edit e : sliceEdits) {
            Rectangle b = e.roi.getBounds();
            ImageProcessor roiMask = e.roi.getMask();
            byte value = e.add ? FG : 0;
            for (int y = Math.max(0, b.y); y < Math.min(h, b.y + b.height); y++) {
                for (int x = Math.max(0, b.x); x < Math.min(w, b.x + b.width); x++) {
                    if (roiMask == null || roiMask.get(x - b.x, y - b.y) != 0) m[y * w + x] = value;
                }
            }
        }
    }

    /**
     * Removes separate objects that look like floating round cells: nearly circular,
     * compact (fills its ellipse) and brighter than the image on average.
     */
    static void excludeRoundBright(byte[] m, float[] raw, int w, int h) {
        int n = w * h;
        double imageMean = mean(raw);
        int[] label = new int[n];
        int[] stack = new int[n];
        int[] region = new int[n];
        int next = 0;
        for (int start = 0; start < n; start++) {
            if (label[start] != 0 || m[start] != FG) continue;
            next++;
            int sp = 0, rn = 0;
            stack[sp++] = start;
            label[start] = next;
            double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sb = 0;
            while (sp > 0) {
                int p = stack[--sp];
                region[rn++] = p;
                int x = p % w, y = p / w;
                sx += x; sy += y; sxx += (double) x * x; syy += (double) y * y; sxy += (double) x * y; sb += raw[p];
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        int nx = x + dx, ny = y + dy;
                        if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
                        int q = ny * w + nx;
                        if (label[q] == 0 && m[q] == FG) { label[q] = next; stack[sp++] = q; }
                    }
                }
            }
            if (rn < 30) continue;
            double mx = sx / rn, my = sy / rn;
            double vxx = sxx / rn - mx * mx, vyy = syy / rn - my * my, vxy = sxy / rn - mx * my;
            double tr = vxx + vyy, det = vxx * vyy - vxy * vxy;
            double disc = Math.sqrt(Math.max(0, tr * tr / 4 - det));
            double l1 = tr / 2 + disc, l2 = Math.max(1e-9, tr / 2 - disc);
            double aspect = Math.sqrt(l1 / l2);
            double fill = rn / (4 * Math.PI * Math.sqrt(l1 * l2));
            double brightness = sb / rn;
            if (aspect < 1.35 && fill > 0.85 && brightness > 1.1 * imageMean) {
                for (int i = 0; i < rn; i++) m[region[i]] = 0;
            }
        }
    }

    /**
     * Sets connected regions of value {@code target} smaller than {@code limit} pixels to
     * {@code replacement}. With {@code keepBorderRegions}, regions touching the image edge are left alone.
     */
    static void flipSmallRegions(byte[] m, int w, int h, byte target, byte replacement, int limit,
                                 boolean eightConnected, boolean keepBorderRegions) {
        int n = w * h;
        boolean[] seen = new boolean[n];
        int[] stack = new int[n];
        int[] region = new int[n];
        for (int start = 0; start < n; start++) {
            if (seen[start] || m[start] != target) continue;
            int sp = 0, rn = 0;
            boolean border = false;
            stack[sp++] = start;
            seen[start] = true;
            while (sp > 0) {
                int p = stack[--sp];
                region[rn++] = p;
                int x = p % w, y = p / w;
                if (x == 0 || y == 0 || x == w - 1 || y == h - 1) border = true;
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        if (dx == 0 && dy == 0) continue;
                        if (!eightConnected && dx != 0 && dy != 0) continue;
                        int nx = x + dx, ny = y + dy;
                        if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
                        int q = ny * w + nx;
                        if (!seen[q] && m[q] == target) {
                            seen[q] = true;
                            stack[sp++] = q;
                        }
                    }
                }
            }
            if (rn < limit && !(keepBorderRegions && border)) {
                for (int i = 0; i < rn; i++) m[region[i]] = replacement;
            }
        }
    }

    /**
     * Removes the bright phase-contrast halo by eroding along the local gradient direction,
     * for at most {@code maxRounds} pixel layers.
     */
    static void haloCorrection(byte[] m, byte[] dir, int w, int h, int maxRounds, BooleanSupplier stale) {
        IntList current = new IntList();
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int p = y * w + x;
                if (m[p] != FG) continue;
                boolean outline = x == 0 || y == 0 || x == w - 1 || y == h - 1
                        || m[p - 1] != FG || m[p + 1] != FG || m[p - w] != FG || m[p + w] != FG;
                if (outline) current.add(p);
            }
        }
        boolean[] visited = new boolean[w * h];
        IntList remove = new IntList();
        int round = 0;
        while (current.size > 0 && round < maxRounds) {
            round++;
            if (stale.getAsBoolean()) throw new Cancelled();
            IntList next = new IntList();
            remove.size = 0;
            for (int i = 0; i < current.size; i++) {
                int p = current.data[i];
                if (visited[p]) continue;
                int x = p % w, y = p / w;
                if (x <= 0 || y <= 0 || x >= w - 1 || y >= h - 1) continue;
                visited[p] = true;
                boolean valid = false;
                int[] cone = CONES[dir[p]];
                for (int k = 0; k < 3; k++) {
                    int[] off = OFFSETS[cone[k] - 1];
                    int q = (y + off[1]) * w + (x + off[0]);
                    if (m[q] == FG) {
                        valid = true;
                        next.add(q);
                    }
                }
                if (valid) remove.add(p);
            }
            for (int i = 0; i < remove.size; i++) m[remove.data[i]] = 0;
            current = next;
        }
    }

    /** Percentage of the analysis area (or whole image when {@code area} is null) covered by cells. */
    static double confluency(ByteProcessor mask, byte[] area) {
        byte[] m = (byte[]) mask.getPixels();
        long on = count(m);
        long total = area == null ? m.length : count(area);
        return total == 0 ? 0 : 100.0 * on / total;
    }

    static long count(byte[] m) {
        long c = 0;
        for (byte b : m) if (b != 0) c++;
        return c;
    }

    static final class IntList {
        int[] data = new int[1024];
        int size;
        void add(int v) {
            if (size == data.length) data = Arrays.copyOf(data, size * 2);
            data[size++] = v;
        }
    }

    // ---------------------------------------------------------------- settings

    static double clamp(double v, double lo, double hi) {
        return Math.max(lo, Math.min(hi, v));
    }
}
