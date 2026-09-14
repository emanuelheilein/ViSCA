# ViSCA

**Visual cell selection and area measurement for immunofluorescence image analysis**

Created by Emanuel Heilein  
First published: 06/2026  
Copyright (C) 2026 Emanuel Heilein

Legal code license: GNU General Public License v3.0  
Legal software license: GNU General Public License v3.0

## Overview

ViSCA is an open-source Python tool for morphology-filtered single-cell area measurement of cultured neonatal rat ventricular cardiomyocytes (NRVCM) in immunofluorescence images. ViSCA processes paired cellbody and nucleus images and combines intensity thresholding, size filtering, mononuclearity assessment, nuclear circularity filtering, nucleus-to-cellbody area ratio filtering, and selectable watershed handling.

The software provides a graphical user interface with interactive preview windows and exports mask panels, CSV files with single-cell measurements, aggregated output tables, and parameter metadata.

## Download ViSCA

Precompiled versions of ViSCA are provided in the **Releases** section of this GitHub repository.

Planned standalone versions are available for:

- Windows 64-bit
- macOS Intel
- macOS Apple Silicon

Alternatively, ViSCA can be run directly from the Python source code as described below.

## Tested environments

ViSCA v1.0 was tested with:

- Python 3.14.7
- NumPy 2.5.3
- SciPy 1.18.1
- pandas 3.0.5
- Matplotlib 3.11.2
- scikit-image 0.26.0
- Pillow 12.3.0
- Tkinter / Tk 9.0

Tested operating environments:

- Windows 11 Pro 64-bit, version 25H2, build 26200.9168, Intel Core i5-13400
- macOS Sequoia 15.1, Apple Silicon M4 Pro
- macOS Catalina 10.15.8, Intel Core i5

Note: Please update this section after final testing of all standalone builds.

## Running ViSCA from source code

### Requirements

The required Python packages are listed in `requirements.txt`.

For the source-code version, keep the following files in the same directory:

```text
ViSCA.py
requirements.txt
logo_ViSCA.png
