# ViSCA

**Visual cell selection and area measurement for immunofluorescence image analysis**

Created by Emanuel Heilein  
First published: 06/2026  
Copyright (C) 2026 Emanuel Heilein

License: GNU General Public License v3.0

---

## Overview

ViSCA is an open-source tool for morphology-filtered single-cell area measurement of cultured neonatal rat ventricular cardiomyocytes (NRVCM) in immunofluorescence images.

ViSCA processes paired cellbody and nucleus images and combines intensity thresholding, size filtering, mononuclearity assessment, nuclear circularity filtering, nucleus-to-cellbody area ratio filtering, and selectable watershed handling.

The software provides a graphical user interface with interactive preview windows and exports mask panels, CSV files with single-cell measurements, aggregated output tables, parameter metadata, and quality-control output files.

---

## Intended use

ViSCA is intended for research use in fluorescence image analysis of cultured cardiomyocytes.

The software was developed for area measurement of cultured neonatal rat ventricular cardiomyocytes stained for cellbody and nuclei. Other cell types, staining protocols, or imaging settings may require parameter adjustment and independent validation.

ViSCA is not intended for clinical diagnosis or medical decision-making.

---

## Download ViSCA without Python installation

For most users, the easiest way to use ViSCA is to download the precompiled app from the latest release:

👉 **[Download ViSCA v1.0](../../releases/latest)**

Choose the file matching your operating system:

| Operating system | Download file |
|---|---|
| Windows 64-bit | `ViSCA_v1.0_Windows_64bit.zip` |
| macOS Intel | `ViSCA_v1.0_macOS_Intel.zip` |
| macOS Apple Silicon (M1/M2/M3/M4) | `ViSCA_v1.0_macOS_AppleSilicon.zip` |

After download, extract the ZIP archive and start:

- `ViSCA.exe` on Windows
- `ViSCA.app` on macOS

No Python installation is required for these precompiled versions.

### Example data

An optional example dataset is available in the release assets:

| File | Description |
|---|---|
| `ViSCA_v1.0_example_data.zip` | Paired example fluorescence images for testing ViSCA |

After extracting the example data ZIP, select:

- `cellbody-channel` as the cellbody image folder
- `nucleus-channel` as the nucleus image folder

The example dataset contains 16 paired fields of view and can be used to test the precompiled app or the source-code version of ViSCA.

### Security note

ViSCA is currently not code-signed. Depending on your operating system settings, Windows or macOS may show a security warning when opening the app for the first time.

On macOS, right-click or Control-click `ViSCA.app`, choose **Open**, and confirm that you want to open the app.

On Windows, a SmartScreen warning may appear. Choose **More info** and then **Run anyway** if you trust the downloaded file.

---

## Quick start

1. Download the ZIP file matching your operating system. Optionally, also download `ViSCA_v1.0_example_data.zip` to test ViSCA with example images.
2. Extract the ZIP archive.
3. Start `ViSCA.exe` on Windows or `ViSCA.app` on macOS.
4. Select the cellbody image folder and the corresponding nucleus image folder.
5. Adjust parameters or use the default settings.
6. Use the preview windows to inspect segmentation and filtering.
7. Start batch processing.
8. Review the generated output folders and CSV files.

---

## Input data

ViSCA requires paired fluorescence microscopy images:

- one image channel showing the cardiomyocyte cellbody
- one image channel showing nuclei

The cellbody and nucleus images should correspond to the same fields of view and should be stored in separate folders.

Example structure:

```text
Experiment_folder/
├── cell/
│   ├── image_01.tif
│   ├── image_02.tif
│   └── image_03.tif
└── nucleus/
    ├── image_01.tif
    ├── image_02.tif
    └── image_03.tif
```

ViSCA pairs cellbody and nuclear images by lexicographically sorted filename order, not by matching filename identifiers. Both folders must therefore contain corresponding images in identical sort order. For example, `image_10.tif` may sort before `image_2.tif`; zero-padded filenames such as `image_02.tif` and `image_10.tif` are recommended. Users must verify channel correspondence before analysis, as incorrect pairing can lead to invalid measurements.

Supported input formats are TIFF (`.tif`, `.tiff`), PNG (`.png`), and JPEG (`.jpg`, `.jpeg`). Lossless grayscale TIFF or PNG images are recommended for quantitative analysis. JPEG compression and RGB-to-grayscale conversion may affect intensity-based segmentation. Because mean-intensity thresholding depends on the image intensity distribution, segmentation thresholds should be checked in the preview windows, particularly when staining intensity, cell density, or image acquisition settings differ.

---

## Output files

For each analyzed image set, ViSCA exports:

- mask panels showing the analysis workflow
- CSV files with single-cell area measurements
- aggregated CSV files
- parameter metadata
- quality-control output files

Typical output includes:

```text
Area_csv/
├── image_01_cell_areas.csv
├── image_02_cell_areas.csv
└── aggregated_cell_areas.csv
```

The area CSV files contain measured cell areas and corresponding nucleus areas for valid selected cells. CSV files use semicolon delimiters (`;`) and decimal commas (`,`). These formatting conventions reflect the native ViSCA v1.0 output. Users importing the files into statistical software should specify the corresponding delimiter and decimal separator.

The exported parameter metadata document the settings used for an analysis run. In ViSCA v1.0, these files are intended for documentation and are not automatically reloaded to replay an analysis. Reproducing a run requires using the same input images, software version, and analysis settings.

---

## Default analysis parameters

The standard ViSCA v1.0 workflow uses the following default parameters:

| Parameter | Default value |
|---|---:|
| Pixel scale | `3.0769 px/µm` |
| Cellbody intensity threshold | `mean` |
| Minimum cellbody area | `100 µm²` |
| Maximum cellbody area | `10,000 µm²` |
| Nucleus intensity threshold | `otsu` |
| Minimum nucleus area | `20 µm²` |
| Maximum nucleus area | `600 µm²` |
| Nuclear circularity threshold | `0.73` |
| Nucleus-to-cellbody area ratio threshold | `0.3` |
| Watershed mode | `Split` |
| Watershed sigma | `6.0` |
| Watershed minimum distance | `30` |
| Watershed h value | `3.0` |

ViSCA retains candidate cell-body objects containing exactly one segmented nucleus fully enclosed by the cell-body mask. Cells overlapping a nuclear object with circularity greater than 0 and less than or equal to the configured circularity threshold (default: 0.73) are excluded. Nuclear circularity is calculated as 4πA/P² using the scikit-image `regionprops.perimeter` estimator. The nucleus-to-cell-body area ratio filter retains cells with a ratio strictly below the configured threshold (default: 0.30); objects with ratios ≥ 0.30 are excluded at the default setting. Minimum and maximum area thresholds are inclusive, and objects touching the image border are excluded.

The default circularity cutoff was selected by visual inspection of NRVCM images. Nuclear circularity and other segmentation results depend on image resolution, staining quality, and morphology. Users should inspect representative previews and adjust parameters when imaging conditions differ.

---

## Watershed modes

ViSCA provides three watershed handling modes:

| Mode | Description |
|---|---|
| `No Watershed` | Cellbody objects are analyzed without watershed-based separation. |
| `Split` | Watershed-separated segments are retained only if they pass downstream size, border, and morphology filters. |
| `Exclude` | Objects affected by watershed separation are excluded from final measurement. |

`Split` is the default mode in ViSCA v1.0.

In Split and Exclude modes, watershed processing is applied to thresholded cell-body objects before cell-size and border filtering. The implementation uses a distance transform, Gaussian smoothing (default σ = 6 pixels), h-maxima processing (default h = 3), and marker detection with a minimum distance of 30 pixels before marker-controlled watershed segmentation. These parameters are adjustable but pixel-dependent and may require adaptation when image resolution or magnification changes.

---

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

- Windows 11 Pro 64-bit, Intel Core i5-13400
- macOS Sequoia 15.1, Apple Silicon M4 Pro
- macOS Catalina 10.15.8, Intel Core i5

---

## Running ViSCA from source code

Users who prefer to run ViSCA directly from the Python source code can do so using the instructions below.

### Required files

The canonical source files are available together in `src/` (`ViSCA.py` and `logo_ViSCA.png`). The repository root also retains the same files for compatibility with the commands below. To run from `src/`, start `python src/ViSCA.py` (Windows) or `python3 src/ViSCA.py` (macOS) from the repository root.

For the root-level launch commands below, keep these files in the same directory:

```text
ViSCA.py
requirements.txt
logo_ViSCA.png
```

### Python requirements

The required Python packages are listed in `requirements.txt`.

Example `requirements.txt` for ViSCA v1.0:

```text
numpy==2.5.3
scipy==1.18.1
pandas==3.0.5
matplotlib==3.11.2
scikit-image==0.26.0
Pillow==12.3.0
```

Tkinter is not listed in `requirements.txt` because it is provided through the Python or system installation.

---

## Windows source-code installation

1. Install Python 3.
2. Download or clone this repository.
3. Open PowerShell in the ViSCA folder.
4. Create and activate a virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

5. Install the required packages:

```powershell
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
```

6. Start ViSCA:

```powershell
python ViSCA.py
```

---

## macOS source-code installation

1. Install Python 3.
2. Download or clone this repository.
3. Open Terminal in the ViSCA folder.
4. Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

5. Install the required packages:

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
```

6. Start ViSCA:

```bash
python ViSCA.py
```

---

## Citation

If you use ViSCA in scientific work, please cite the corresponding publication once available.

Citation metadata are provided in:

```text
CITATION.cff
```

---

## License

ViSCA is distributed under the GNU General Public License v3.0.

See `LICENSE.txt` for the complete GPL v3 license text (`LICENSE` is retained for GitHub license detection).

---

## Contact

For questions, feedback, or bug reports, please contact:

Emanuel Heilein  
emanuel@heilein.de
