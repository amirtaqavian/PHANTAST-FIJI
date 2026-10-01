# PHANTAST Live — interactive confluency measurement for Fiji

**PHANTAST Live** is a Fiji/ImageJ plugin that measures **cell confluency in phase-contrast microscopy images**, with a live preview, image adjustments, analysis-area selection and manual corrections.

It is based on the original [PHANTAST for FIJI](https://github.com/nicjac/PHANTAST-FIJI) by Nicolas Jaccard (UCL) and uses the same segmentation method (local contrast thresholding with halo correction, Jaccard *et al.*, 2014), re-implemented with an interactive interface. The original PHANTAST plugin code is kept unchanged in this repository.

## Download and install

1. Download **`PHANTAST_Live-1.0.jar`** from the [latest release](https://github.com/a-LittleMonk/PHANTAST-FIJI/releases/latest).
2. Copy it into the **`plugins`** folder of your Fiji installation (or drag it onto the Fiji toolbar and save it in `plugins`).
3. Restart Fiji.
4. Open an image and run **Plugins › Segmentation › PHANTAST Live**.

Requires **Fiji** (or ImageJ 1.54 or newer). It works on 8-bit, 16-bit, 32-bit and RGB images, single images and stacks. PHANTAST Live can be installed next to the original PHANTAST plugin; they do not interfere.

## How to use it

The settings window has a preset bar at the top, five tabs in workflow order, and the live **confluency** at the bottom.

| Tab | What you do there |
| --- | --- |
| **1. Area** | Choose what to analyse: the **whole image**, a **selection drawn on the image** (press *Set area from selection*), or **auto-detect the circular field of view** of photos taken through the eyepiece (with an edge margin). |
| **2. Image adjustments** | Brightness, contrast, clarity, gamma, local contrast (CLAHE), median and Gaussian smoothing, rolling-ball background correction and even illumination. Changes show **instantly**; press **Apply to detection** to detect cells on the adjusted image, or **Default** to reset. Cell detection is hidden on this tab by default so you can see the cells clearly. |
| **3. Cell detection** | Live preview options (yellow outline, green fill, live black & white mask window), then **sigma** and **epsilon** (the PHANTAST parameters), halo correction strength, minimum cell size, hole filling, grow/shrink, and *exclude round bright cells*. |
| **4. Manual correction** | Draw a selection on the image and press **Add selection to cells** or **Remove selection**; undo or clear edits. |
| **5. Output** | What happens when you click **OK**: add a row to the Results table, draw the yellow cell outline on the image, create a black & white mask, process all slices of a stack. |

**Presets:** the plugin always opens on **Default**. Change any settings and press **Save preset** to store them under a name. When a saved preset is selected and modified, *Save preset* asks whether to **update the preset**, **save a new preset** or **cancel**. Presets are stored in the ImageJ preferences folder (`PHANTAST_Live_presets.properties`).

**Results table:** each row records the image, confluency (%), analysis area and its size, and every setting used (sigma, epsilon, halo correction, size and hole limits, grow/shrink, round-cell exclusion, number of manual edits and any image adjustments applied to detection), so results can be reported and reproduced.

**Batch / macro use:** run the plugin with a saved preset without opening the window:

```
run("PHANTAST Live", "preset=[My preset]");
```

### Tips

* Use the **same settings (preset) for every image you compare**, including controls.
* Image adjustments applied to detection can change the result substantially — check the yellow outline after pressing *Apply*. Strong clarity or CLAHE make background texture look like cells; raise epsilon if that happens.
* *Exclude round bright cells* removes separate, round floating cells; round cells touching other cells are not removed — use *Remove selection* for those.

## Building from source

The plugin is a single Java file, [`PHANTAST_Live/src/PHANTAST_Live.java`](PHANTAST_Live/src/PHANTAST_Live.java), that depends only on the core ImageJ library. To build it with the Java compiler bundled with Fiji:

* **Windows:** `PHANTAST_Live\build.bat "C:\path\to\Fiji"`
* **macOS / Linux** (needs `javac` 8+ on the PATH): `./PHANTAST_Live/build.sh /path/to/Fiji.app`

This creates `PHANTAST_Live/PHANTAST_Live-1.0.jar`. The source code is also included inside the released jar.

## Repository layout

| Path | Contents |
| --- | --- |
| `PHANTAST_Live/` | PHANTAST Live source code, menu configuration and build scripts |
| `src/`, `pom.xml` | Original PHANTAST for FIJI plugin (unchanged) |
| `LICENSE` | License (BSD 3-clause, from the original PHANTAST) |

## Citation

If you use this plugin in your research, please cite the original PHANTAST publication:

> Jaccard N, Griffin LD, Keser A, Macown RJ, Super A, Veraitch FS, Szita N. **Automated method for the rapid and precise estimation of adherent cell culture characteristics from phase contrast microscopy images.** *Biotechnology and Bioengineering* 2014; 111(3): 504–517. doi:[10.1002/bit.25115](https://doi.org/10.1002/bit.25115)

and mention that confluency was measured with PHANTAST Live together with the settings or preset used.

## License and credits

PHANTAST Live is a modified version of PHANTAST for FIJI and is distributed under the same BSD 3-clause license (see [`LICENSE`](LICENSE)); the original copyright notice is retained in the source code.

* Original PHANTAST and PHANTAST for FIJI: Nicolas Jaccard, Department of Biochemical Engineering, UCL — <https://github.com/nicjac/PHANTAST-FIJI>
* 2017 FIJI plugin improvements: Olivier Burri, BIOP, EPFL
* PHANTAST Live (interactive interface, image adjustments, analysis areas, presets, manual corrections): a-LittleMonk
