# PHANTAST_Live — Enhanced PHANTAST for FIJI

**PHANTAST_Live** is a modified version of the original **PHANTAST for FIJI** plugin for phase-contrast microscopy image analysis.

This version builds on the original PHANTAST segmentation workflow and introduces additional features focused on **live detection, image adjustment, and easier Region of Interest (ROI) selection**.

## Getting Started

Download the latest version of **PHANTAST_Live** from the release page:

**[Download PHANTAST_Live](https://github.com/amirtaqavian/PHANTAST-FIJI/releases)**

After downloading, install the plugin in **Fiji/ImageJ** and start using PHANTAST_Live for phase-contrast microscopy image analysis.

## What's New

### Live Detection

* Added **live detection** functionality for interactive image analysis.
* Detection results can be updated while adjusting the image and analysis settings.
* Provides a more immediate way to evaluate segmentation results.

### Image Adjustment

* Added image adjustment functionality to make preprocessing easier.
* Users can modify the image before running detection to obtain more suitable segmentation results.
* Provides greater control when working with images with different contrast, brightness, or background characteristics.

### Easier ROI Selection

* Improved the workflow for selecting the **Region of Interest (ROI)**.
* Users can more easily define and adjust the area that should be analyzed.
* This helps exclude unwanted areas and focus detection on the relevant part of the microscopy image.

## Why PHANTAST_Live?

The original PHANTAST was developed to provide automated segmentation and quantitative analysis of phase-contrast microscopy images, including cell confluency and cell density estimation.

**PHANTAST_Live extends this workflow by making the analysis more interactive**, allowing users to adjust the input image and ROI while directly observing the resulting detection.

## Main Changes

| Feature                           | Original PHANTAST | PHANTAST_Live |
| --------------------------------- | ----------------- | ------------- |
| Phase-contrast image segmentation | ✓                 | ✓             |
| Cell confluency analysis          | ✓                 | ✓             |
| Cell density estimation           | ✓                 | ✓             |
| Live detection                    | —                 | ✓             |
| Interactive image adjustment      | Limited           | ✓             |
| Easier ROI selection              | Limited           | ✓             |
| Interactive analysis workflow     | —                 | ✓             |


## Note

PHANTAST_Live is an independent modification of the original PHANTAST-FIJI plugin. It is intended to provide additional functionality while maintaining the core purpose of PHANTAST: quantitative analysis of phase-contrast microscopy images.
