# ViSCA
**Visual cell selection and area measurement for immunofluorescence image analysis**

Created by Emanuel Heilein
First published: 06/2026
Copyright (C) 2026 Emanuel Heilein

Legal Code License: GNU General Public License v3.0  
Legal Software License: GNU General Public License v3.0

## Download ViSCA

Precompiled versions of ViSCA are available for:

- Windows 64-bit
- macOS Intel
- macOS Apple Silicon

The precompiled versions are provided in the **Releases** section of this GitHub repository.
Alternatively, ViSCA can be run directly from the Python source code as described below.

## Running ViSCA from source code

### Requirements

ViSCA v1.0 was tested with Python 3.14.7 and the package versions specified in `requirements.txt`.

Main dependencies:

- NumPy
- SciPy
- pandas
- Matplotlib
- scikit-image
- Pillow
- Tkinter

The required Python packages are listed in `requirements.txt`.

For the source-code version, keep the following files in the same directory:

```text
visca.py
requirements.txt
logo_ViSCA.png
```

### Windows

1. Install a current version of Python 3 from:
   https://www.python.org/downloads/windows/
   During installation, enable:

   `Add python.exe to PATH`

2. Download the ViSCA source files and place `visca.py`, `requirements.txt`, and `logo_ViSCA.png` in the same directory.

3. Open **PowerShell** or **Command Prompt** in this directory.

4. Create a virtual Python environment:

```powershell
py -m venv venv
```

5. Install the required packages:

```powershell
venv\Scripts\python.exe -m pip install --upgrade pip
venv\Scripts\python.exe -m pip install -r requirements.txt
```

6. Start ViSCA:

```powershell
venv\Scripts\python.exe visca.py
```

### macOS

1. Install a current version of Python 3 from:
   https://www.python.org/downloads/macos/

2. Download the ViSCA source files and place `visca.py`, `requirements.txt`, and `logo_ViSCA.png` in the same directory.

3. Open **Terminal** in this directory.

4. Create a virtual Python environment:

```bash
python3 -m venv venv
```

5. Install the required packages:

```bash
./venv/bin/python -m pip install --upgrade pip
./venv/bin/python -m pip install -r requirements.txt
```

6. Start ViSCA:

```bash
./venv/bin/python visca.py
```

## License

ViSCA is free and open-source software distributed under the **GNU General Public License v3.0**.
See the `LICENSE` file in this repository for the complete license text.
## Citation

If you use ViSCA in scientific work, please cite the corresponding publication.
The complete citation and DOI will be added after publication.
