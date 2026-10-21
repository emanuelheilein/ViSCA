# ViSCA - Visual cell selection and area measurement
# Version 1.0
# Copyright (C) 2026 Emanuel Heilein
#
# SPDX-License-Identifier: GPL-3.0-only
#
# This program is free software licensed under the
# GNU General Public License version 3.

import os
import sys
import math
import multiprocessing
import threading
import time
import tkinter as tk
from math import pi
from tkinter import ttk, messagebox, filedialog
from PIL import Image, ImageTk, ImageOps

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.dirname(os.path.abspath(__file__))

    return os.path.join(base_path, relative_path)

def load_analysis_modules():
    global plt, np, pd
    global LinearSegmentedColormap
    global ndimage
    global filters, io, measure
    global peak_local_max, gaussian, h_maxima, watershed, find_boundaries

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    import pandas as pd
    from scipy import ndimage
    from skimage import filters, io, measure
    from skimage.feature import peak_local_max
    from skimage.filters import gaussian
    from skimage.morphology import h_maxima
    from skimage.segmentation import watershed, find_boundaries

def show_loading_screen_and_load():

    from PIL import ImageEnhance

    splash = tk.Tk()
    splash.title("ViSCA")
    splash.resizable(False, False)
    splash.configure(bg="white")

    window_width = 800
    window_height = 720

    screen_width = splash.winfo_screenwidth()
    screen_height = splash.winfo_screenheight()

    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    splash.geometry(f"{window_width}x{window_height}+{x}+{y}")

    frame = tk.Frame(splash, bg="white")
    frame.pack(fill="both", expand=True)

    try:
        logo_image = Image.open(resource_path("logo_ViSCA.png"))
        logo_image = ImageOps.exif_transpose(logo_image)
        logo_image = logo_image.convert("RGB")
        logo_image = ImageEnhance.Color(logo_image).enhance(2)

        max_width = 600
        max_height = 400

        logo_image.thumbnail(
            (max_width, max_height),
            Image.Resampling.LANCZOS
        )

        logo = ImageTk.PhotoImage(logo_image)

        logo_label = tk.Label(
            frame,
            image=logo,
            bg="white",
            borderwidth=0
        )
        logo_label.image = logo
        logo_label.pack(pady=(35, 25))

    except Exception as exc:
        logo_label = tk.Label(
            frame,
            text="ViSCA",
            bg="white",
            fg="#4FD354",
            font=("Arial", 48, "bold")
        )
        logo_label.pack(pady=(120, 60))
        print("Logo could not be loaded:", exc)

    info_label = tk.Label(
        frame,
        text=(
            "ViSCA - Visual cell selection and area measurement\n"
            "Version 1.0\n"
            "Developed by Emanuel Heilein\n"
            "Copyright (C) 2026 Emanuel Heilein\n"
            "Licensed under the GNU General Public License v3.0"
        ),
        bg="white",
        fg="black",
        justify="center",
        font=("Arial", 14)
    )
    info_label.pack(pady=(0, 35))

    progress = ttk.Progressbar(
        frame,
        mode="determinate",
        length=420,
        maximum=100,
        value=0
    )
    progress.pack()

    status_label = tk.Label(
        frame,
        text="Loading ViSCA... 0%",
        bg="white",
        fg="black",
        font=("Arial", 11)
    )
    status_label.pack(pady=(12, 0))

    splash.update_idletasks()
    splash.update()

    start_time = time.time()
    estimated_seconds = 30
    loading_done = {"value": False}

    def load_in_background():
        load_analysis_modules()
        loading_done["value"] = True

    def update_progress_bar():
        elapsed = time.time() - start_time
        percent = min(int((elapsed / estimated_seconds) * 100), 99)

        if loading_done["value"]:
            progress["value"] = 100
            status_label.config(text="Loading ViSCA... 100%")
            splash.update_idletasks()
            splash.after(300, splash.destroy)
            return

        progress["value"] = percent
        status_label.config(text=f"Loading ViSCA... {percent}%")
        splash.after(100, update_progress_bar)

    threading.Thread(target=load_in_background, daemon=True).start()

    update_progress_bar()

    splash.mainloop()

DEFAULT_PARAMETERS = {
    "scale": 3.0769,
    "intensity_threshold_cell": "mean",
    "min_size_cell": 100,
    "max_size_cell": 10000,
    "intensity_threshold_nucleus": "otsu",
    "min_size_nucleus": 20,
    "max_size_nucleus": 600,
    "circularity_threshold": 0.73,
    "ratio_threshold": 0.3,
    "watershed_mode": "Split",
    "watershed_sigma": 6.0,
    "watershed_min_distance": 30,
    "watershed_h_value": 3.0,
}

IMAGE_EXTENSIONS = (".tif", ".tiff", ".png", ".jpg", ".jpeg")
NUM_PROCESSES = 4


# Helper functions

def stepped_slider_sequence(x):
    return 10 ** (x // 9) * (x % 9 + 1)


def slider_position_value(x):
    n_z = math.log10(x) // 1
    r = n_z * 9 + (x / 10 ** n_z) - 1
    return round(r)

def ensure_directory(path):
    os.makedirs(path, exist_ok=True)


def list_image_files(folders):
    files = []
    for folder in folders:
        for filename in os.listdir(folder):
            if filename.lower().endswith(IMAGE_EXTENSIONS) and not filename.startswith("."):
                files.append(os.path.join(folder, filename))
    files.sort()
    return files

def validate_parameters(parameters):
    if parameters["min_size_cell"] >= parameters["max_size_cell"]:
        raise ValueError("min_size_cell must be less than max_size_cell")

    if parameters["min_size_nucleus"] >= parameters["max_size_nucleus"]:
        raise ValueError("min_size_nucleus must be less than max_size_nucleus")

    if not (0 < parameters["circularity_threshold"] <= 1):
        raise ValueError("circularity must be between 0 and 1")

    if parameters["ratio_threshold"] < 0:
        raise ValueError("ratio_threshold must not be negative")

    if parameters["scale"] <= 0:
        raise ValueError("scale must be greater than 0")

    if parameters.get("watershed_mode") not in {"No Watershed", "Split", "Exclude"}:
        raise ValueError("watershed_mode must be 'No Watershed', 'Split', or 'Exclude'")

    if parameters["watershed_sigma"] < 0:
        raise ValueError("watershed_sigma must not be negative")

    if parameters["watershed_min_distance"] < 1:
        raise ValueError("watershed_min_distance must be at least 1")

    if parameters["watershed_h_value"] < 0:
        raise ValueError("watershed_h_value must not be negative")


def add_placeholder(entry, placeholder):
    entry_color = entry.cget("foreground")
    entry._default_foreground = entry_color

    def on_focus_in(_event, placeholder_text=placeholder):
        if entry.get() == placeholder_text and entry.cget("foreground") == "grey":
            entry.delete(0, "end")
            entry.config(foreground=entry._default_foreground)

    def on_focus_out(_event, placeholder_text=placeholder):
        if not entry.get():
            entry.insert(0, placeholder_text)
            entry.config(foreground="grey")
        else:
            entry.config(foreground=entry._default_foreground)

    entry.insert(0, placeholder)
    entry.config(foreground="grey")
    entry.bind("<FocusIn>", on_focus_in)
    entry.bind("<FocusOut>", on_focus_out)


def set_entry_value(entry, value):
    entry.delete(0, tk.END)
    entry.insert(0, str(value))
    default_fg = getattr(entry, "_default_foreground", entry.cget("foreground"))
    entry.config(foreground=default_fg)


def enable_arrow_key_step(scale_widget):
    scale_widget.configure(takefocus=True)

    def focus_slider(_event):
        scale_widget.focus_set()

    def step_left(_event):
        current = float(scale_widget.get())
        step = float(scale_widget.cget("resolution"))
        minimum = float(scale_widget.cget("from"))
        scale_widget.set(max(minimum, current - step))
        return "break"

    def step_right(_event):
        current = float(scale_widget.get())
        step = float(scale_widget.cget("resolution"))
        maximum = float(scale_widget.cget("to"))
        scale_widget.set(min(maximum, current + step))
        return "break"

    scale_widget.bind("<Button-1>", focus_slider)
    scale_widget.bind("<Left>", step_left)
    scale_widget.bind("<Right>", step_right)

required_color = "#4FD354"

# Image analysis functions

def intensity_threshold(image_path, min_size, max_size, intensity_threshold_value):
    image = io.imread(image_path)
    image_gray = image.mean(axis=2) if image.ndim == 3 else image

    if intensity_threshold_value == "otsu":
        intensity_threshold_value = filters.threshold_otsu(image_gray)
    elif intensity_threshold_value == "mean":
        intensity_threshold_value = image_gray.mean()

    try:
        intensity_threshold_value = float(intensity_threshold_value)
    except ValueError as exc:
        raise ValueError("Invalid intensity threshold value") from exc

    binary_mask = image_gray > intensity_threshold_value
    binary_mask_filled = ndimage.binary_fill_holes(binary_mask)
    labeled_mask = measure.label(binary_mask_filled)
    regions = measure.regionprops(labeled_mask)

    final_mask = np.zeros_like(binary_mask_filled, dtype=bool)

    for region in regions:
        if (
            min_size <= region.area <= max_size
            and all(x > 0 for x in region.bbox[:2])
            and all(x < y for x, y in zip(region.bbox[2:], binary_mask_filled.shape))
        ):
            final_mask[region.coords[:, 0], region.coords[:, 1]] = True

    return ~final_mask, intensity_threshold_value

def intensity_threshold_with_watershed_split(
    image_path,
    min_size,
    max_size,
    intensity_threshold_value,
    sigma=6.0,
    min_distance=30,
    h_value=3.0,
):
    image = io.imread(image_path)
    image_gray = image.mean(axis=2) if image.ndim == 3 else image

    if intensity_threshold_value == "otsu":
        intensity_threshold_value = filters.threshold_otsu(image_gray)
    elif intensity_threshold_value == "mean":
        intensity_threshold_value = image_gray.mean()

    try:
        intensity_threshold_value = float(intensity_threshold_value)
    except ValueError as exc:
        raise ValueError("Invalid intensity threshold value") from exc

    binary_mask = image_gray > intensity_threshold_value
    binary_mask_filled = ndimage.binary_fill_holes(binary_mask)

    # Original connected objects before watershed.
    # These are the objects that should remain unchanged if the watershed does not really split them.
    pre_watershed_labels = measure.label(binary_mask_filled)

    original_object_areas_pixels = {
        int(region.label): int(region.area)
        for region in measure.regionprops(pre_watershed_labels)
    }

    pre_watershed_regions = {
        int(region.label): region
        for region in measure.regionprops(pre_watershed_labels)
    }

    # For plotting: thresholded cells before watershed and before size filtering.
    pre_watershed_mask = ~binary_mask_filled

    distance = ndimage.distance_transform_edt(binary_mask_filled)
    smoothed_distance = gaussian(distance, sigma=sigma)
    h_maxima_transformed = h_maxima(smoothed_distance, h=h_value)

    local_maxi = peak_local_max(
        h_maxima_transformed,
        footprint=np.ones((3, 3)),
        min_distance=min_distance,
        labels=binary_mask_filled,
        exclude_border=False,
    )

    markers = np.zeros_like(binary_mask_filled, dtype=int)

    if len(local_maxi) > 0:
        markers[tuple(local_maxi.T)] = np.arange(1, len(local_maxi) + 1)
        watershed_labels = watershed(
            -smoothed_distance,
            markers,
            mask=binary_mask_filled,
        )
    else:
        watershed_labels = measure.label(binary_mask_filled)

    watershed_regions = {
        int(region.label): region
        for region in measure.regionprops(watershed_labels)
    }

    # Count all raw watershed segments per original pre-watershed object,
    # before size and border filtering.
    raw_original_to_segments = {}

    for raw_region in measure.regionprops(watershed_labels):
        raw_rows = raw_region.coords[:, 0]
        raw_cols = raw_region.coords[:, 1]

        raw_original_ids = pre_watershed_labels[raw_rows, raw_cols]
        raw_original_ids = raw_original_ids[raw_original_ids > 0]

        if len(raw_original_ids) > 0:
            raw_original_id = int(np.bincount(raw_original_ids).argmax())
        else:
            raw_original_id = 0

        if raw_original_id > 0:
            raw_original_to_segments.setdefault(raw_original_id, []).append(
                int(raw_region.label)
            )

    # Analysis label mask:
    # - If an original object has only one raw watershed segment, keep the original object mask.
    # - If an original object has multiple raw watershed segments, use the watershed segment masks.
    filtered_watershed_labels = np.zeros_like(watershed_labels, dtype=int)
    watershed_segment_info = {}
    new_label = 1

    for original_object_id, original_region in pre_watershed_regions.items():
        raw_segment_labels = sorted(
            raw_original_to_segments.get(original_object_id, [])
        )

        raw_number_of_segments = len(raw_segment_labels)
        if raw_number_of_segments < 1:
            raw_number_of_segments = 1

        original_area_pixels = original_object_areas_pixels.get(
            original_object_id,
            int(original_region.area),
        )

        # Case 1: Object was not split by watershed.
        # Keep exactly the original pre-watershed object area.
        if raw_number_of_segments == 1:
            if (
                min_size <= original_region.area <= max_size
                and all(x > 0 for x in original_region.bbox[:2])
                and all(
                    x < y
                    for x, y in zip(
                        original_region.bbox[2:],
                        binary_mask_filled.shape,
                    )
                )
            ):
                rows = original_region.coords[:, 0]
                cols = original_region.coords[:, 1]

                filtered_watershed_labels[rows, cols] = new_label

                watershed_segment_info[new_label] = {
                    "original_object_id": original_object_id,
                    "original_object_area_pixels": original_area_pixels,
                    "raw_number_of_segments_in_original_object": raw_number_of_segments,
                }

                new_label += 1

            continue

        # Case 2: Object was split by watershed.
        # Keep watershed segments as separate measurable objects.
        for raw_segment_label in raw_segment_labels:
            region = watershed_regions.get(raw_segment_label)

            if region is None:
                continue

            if (
                min_size <= region.area <= max_size
                and all(x > 0 for x in region.bbox[:2])
                and all(
                    x < y
                    for x, y in zip(
                        region.bbox[2:],
                        binary_mask_filled.shape,
                    )
                )
            ):
                rows = region.coords[:, 0]
                cols = region.coords[:, 1]

                filtered_watershed_labels[rows, cols] = new_label

                watershed_segment_info[new_label] = {
                    "original_object_id": original_object_id,
                    "original_object_area_pixels": original_area_pixels,
                    "raw_number_of_segments_in_original_object": raw_number_of_segments,
                }

                new_label += 1

    # Count how many retained watershed segments came from each original object.
    original_to_segments = {}

    for segment_label, info in watershed_segment_info.items():
        original_id = info["original_object_id"]
        original_to_segments.setdefault(original_id, []).append(segment_label)

    for segment_label, info in watershed_segment_info.items():
        original_id = info["original_object_id"]
        info["number_of_segments_in_original_object"] = len(
            original_to_segments.get(original_id, [])
        )

    # Plot mask only: show internal watershed boundaries as visual lines.
    all_boundaries = find_boundaries(filtered_watershed_labels, mode="inner")
    outer_boundaries = find_boundaries(filtered_watershed_labels > 0, mode="inner")
    internal_boundaries = all_boundaries & ~outer_boundaries

    plot_mask = filtered_watershed_labels > 0
    plot_mask_with_lines = plot_mask.copy()
    plot_mask_with_lines[internal_boundaries] = False

    return (
        ~plot_mask_with_lines,
        intensity_threshold_value,
        pre_watershed_mask,
        filtered_watershed_labels,
        watershed_segment_info,
    )


def threshold_cellbody_without_size_filter(image_path, intensity_threshold_value):
    """
    Cellbody thresholding only.
    Returns the ViSCA-style inverted threshold mask plus the filled binary cellbody mask.
    No size or border filtering is applied here.
    """
    image = io.imread(image_path)
    image_gray = image.mean(axis=2) if image.ndim == 3 else image

    if intensity_threshold_value == "otsu":
        intensity_threshold_value = filters.threshold_otsu(image_gray)
    elif intensity_threshold_value == "mean":
        intensity_threshold_value = image_gray.mean()

    try:
        intensity_threshold_value = float(intensity_threshold_value)
    except ValueError as exc:
        raise ValueError("Invalid intensity threshold value") from exc

    binary_mask = image_gray > intensity_threshold_value
    binary_mask_filled = ndimage.binary_fill_holes(binary_mask)
    pre_watershed_mask = ~binary_mask_filled

    return pre_watershed_mask, intensity_threshold_value, binary_mask_filled


def build_watershed_labels_without_size_filter(
    binary_mask_filled,
    sigma=6.0,
    min_distance=30,
    h_value=3.0,
):
    """
    Apply watershed to thresholded cellbody objects before size filtering.

    If an original thresholded object is not split, its original object mask is kept.
    If it is split, the raw watershed segments are used.
    No size or border filtering is applied in this function.
    """
    pre_watershed_labels = measure.label(binary_mask_filled)

    original_object_areas_pixels = {
        int(region.label): int(region.area)
        for region in measure.regionprops(pre_watershed_labels)
    }

    pre_watershed_regions = {
        int(region.label): region
        for region in measure.regionprops(pre_watershed_labels)
    }

    distance = ndimage.distance_transform_edt(binary_mask_filled)
    smoothed_distance = gaussian(distance, sigma=sigma)
    h_maxima_transformed = h_maxima(smoothed_distance, h=h_value)

    local_maxi = peak_local_max(
        h_maxima_transformed,
        footprint=np.ones((3, 3)),
        min_distance=min_distance,
        labels=binary_mask_filled,
        exclude_border=False,
    )

    markers = np.zeros_like(binary_mask_filled, dtype=int)

    if len(local_maxi) > 0:
        markers[tuple(local_maxi.T)] = np.arange(1, len(local_maxi) + 1)
        watershed_labels = watershed(
            -smoothed_distance,
            markers,
            mask=binary_mask_filled,
        )
    else:
        watershed_labels = measure.label(binary_mask_filled)

    watershed_regions = {
        int(region.label): region
        for region in measure.regionprops(watershed_labels)
    }

    raw_original_to_segments = {}

    for raw_region in measure.regionprops(watershed_labels):
        raw_rows = raw_region.coords[:, 0]
        raw_cols = raw_region.coords[:, 1]

        raw_original_ids = pre_watershed_labels[raw_rows, raw_cols]
        raw_original_ids = raw_original_ids[raw_original_ids > 0]

        if len(raw_original_ids) > 0:
            raw_original_id = int(np.bincount(raw_original_ids).argmax())
        else:
            raw_original_id = 0

        if raw_original_id > 0:
            raw_original_to_segments.setdefault(raw_original_id, []).append(
                int(raw_region.label)
            )

    raw_cellbody_labels = np.zeros_like(watershed_labels, dtype=int)
    watershed_segment_info = {}
    new_label = 1

    for original_object_id, original_region in pre_watershed_regions.items():
        raw_segment_labels = sorted(
            raw_original_to_segments.get(original_object_id, [])
        )

        raw_number_of_segments = len(raw_segment_labels)
        if raw_number_of_segments < 1:
            raw_number_of_segments = 1

        original_area_pixels = original_object_areas_pixels.get(
            original_object_id,
            int(original_region.area),
        )

        if raw_number_of_segments == 1:
            rows = original_region.coords[:, 0]
            cols = original_region.coords[:, 1]
            raw_cellbody_labels[rows, cols] = new_label

            watershed_segment_info[new_label] = {
                "original_object_id": original_object_id,
                "original_object_area_pixels": original_area_pixels,
                "raw_number_of_segments_in_original_object": raw_number_of_segments,
                "number_of_segments_in_original_object": 1,
            }

            new_label += 1
            continue

        for raw_segment_label in raw_segment_labels:
            region = watershed_regions.get(raw_segment_label)
            if region is None:
                continue

            rows = region.coords[:, 0]
            cols = region.coords[:, 1]
            raw_cellbody_labels[rows, cols] = new_label

            watershed_segment_info[new_label] = {
                "original_object_id": original_object_id,
                "original_object_area_pixels": original_area_pixels,
                "raw_number_of_segments_in_original_object": raw_number_of_segments,
                "number_of_segments_in_original_object": raw_number_of_segments,
            }

            new_label += 1

    return raw_cellbody_labels, watershed_segment_info


def size_filter_label_image(label_image, segment_info, min_size, max_size):
    """
    Apply cell-size and border filtering to an already labeled cellbody image.
    Labels are relabeled consecutively after filtering.
    Returns filtered labels, relabeled metadata, and a mask of removed objects/segments.
    """
    label_image = np.asarray(label_image)
    filtered_labels = np.zeros_like(label_image, dtype=int)
    filtered_info = {}
    removed_mask = np.zeros_like(label_image, dtype=bool)
    new_label = 1

    for region in measure.regionprops(label_image):
        old_label = int(region.label)

        keep = (
            min_size <= region.area <= max_size
            and all(x > 0 for x in region.bbox[:2])
            and all(x < y for x, y in zip(region.bbox[2:], label_image.shape))
        )

        rows = region.coords[:, 0]
        cols = region.coords[:, 1]

        if keep:
            filtered_labels[rows, cols] = new_label
            info = dict(segment_info.get(old_label, {}))
            if not info:
                info = {
                    "original_object_id": old_label,
                    "original_object_area_pixels": int(region.area),
                    "raw_number_of_segments_in_original_object": 1,
                    "number_of_segments_in_original_object": 1,
                }
            filtered_info[new_label] = info
            new_label += 1
        else:
            removed_mask[rows, cols] = True

    original_to_segments = {}
    for segment_label, info in filtered_info.items():
        original_id = info.get("original_object_id", segment_label)
        original_to_segments.setdefault(original_id, []).append(segment_label)

    for segment_label, info in filtered_info.items():
        original_id = info.get("original_object_id", segment_label)
        info["number_of_segments_in_original_object"] = len(
            original_to_segments.get(original_id, [])
        )

    return filtered_labels, filtered_info, removed_mask


def subset_label_image_by_predicate(label_image, segment_info, predicate):
    """Keep labels whose metadata fulfill predicate(info), then relabel consecutively."""
    out = np.zeros_like(label_image, dtype=int)
    out_info = {}
    new_label = 1

    for region in measure.regionprops(label_image):
        old_label = int(region.label)
        info = dict(segment_info.get(old_label, {}))
        if not predicate(info):
            continue

        rows = region.coords[:, 0]
        cols = region.coords[:, 1]
        out[rows, cols] = new_label
        out_info[new_label] = info
        new_label += 1

    return out, out_info


def watershed_title_number(value):
    """Compact number formatting for mask titles."""
    try:
        value = float(value)
        if value.is_integer():
            return str(int(value))
        return f"{value:g}"
    except Exception:
        return str(value)


def cell_size_range_title(min_size_cell, max_size_cell):
    return f"Cell size filter ({int(min_size_cell)}-{int(max_size_cell)}µm²)"


def run_cellbody_stage(
    image_path,
    min_size_pixels,
    max_size_pixels,
    intensity_threshold_value,
    watershed_mode="No Watershed",
    sigma=6.0,
    min_distance=30,
    h_value=3.0,
):
    """
    Run the cellbody branch in separated display/analysis stages.

    Display stages:
    1. intensity threshold only
    2. watershed/no-watershed display without size filtering
    3. cell size + border filter

    Analysis stages:
    No Watershed: threshold -> size/border filter
    Split: threshold -> watershed -> size/border filter on segments
    Exclude: threshold -> watershed detection -> remove split-affected original objects -> size/border filter
    """
    pre_watershed_mask, threshold_value, binary_mask_filled = threshold_cellbody_without_size_filter(
        image_path,
        intensity_threshold_value,
    )

    if watershed_mode == "No Watershed":
        raw_labels = measure.label(binary_mask_filled)
        raw_info = {}
        watershed_display = labels_to_gray_display(raw_labels, raw_info)

        size_filtered_labels, size_filtered_info, size_removed_mask = size_filter_label_image(
            raw_labels,
            raw_info,
            min_size_pixels,
            max_size_pixels,
        )
        analysis_labels = size_filtered_labels
        analysis_info = size_filtered_info
        size_stage_removed_mask = size_removed_mask
        size_stage_display = labels_to_gray_display(
            analysis_labels,
            analysis_info,
        )
        final_outline_reference_labels = raw_labels

        return {
            "mask_cell": ~(analysis_labels > 0),
            "threshold_cell": threshold_value,
            "pre_watershed_cell_mask": pre_watershed_mask,
            "raw_labels": raw_labels,
            "raw_info": raw_info,
            "watershed_display": watershed_display,
            "size_filtered_labels": size_filtered_labels,
            "size_filtered_info": size_filtered_info,
            "size_removed_mask": size_removed_mask,
            "analysis_labels": analysis_labels,
            "analysis_info": analysis_info,
            "size_stage_display": size_stage_display,
            "stage_display": size_stage_display,
            "final_outline_reference_labels": final_outline_reference_labels,
        }

    raw_labels, raw_info = build_watershed_labels_without_size_filter(
        binary_mask_filled,
        sigma=sigma,
        min_distance=min_distance,
        h_value=h_value,
    )
    watershed_display = labels_to_gray_display(raw_labels, raw_info)

    size_filtered_labels, size_filtered_info, size_removed_mask = size_filter_label_image(
        raw_labels,
        raw_info,
        min_size_pixels,
        max_size_pixels,
    )

    if watershed_mode == "Split":
        analysis_labels = size_filtered_labels
        analysis_info = size_filtered_info
        size_stage_removed_mask = size_removed_mask
        size_stage_display = labels_to_gray_display(
            analysis_labels,
            analysis_info,
        )
    elif watershed_mode == "Exclude":
        not_split_labels, not_split_info = subset_label_image_by_predicate(
            raw_labels,
            raw_info,
            lambda info: int(info.get("raw_number_of_segments_in_original_object", 1)) <= 1,
        )
        analysis_labels, analysis_info, size_removed_after_exclude = size_filter_label_image(
            not_split_labels,
            not_split_info,
            min_size_pixels,
            max_size_pixels,
        )
        # In Exclude mode, watershed-affected original objects are removed before
        # downstream filters. Objects removed by watershed exclusion, size, or
        # border filtering are intentionally hidden in the size-filter panel.
        size_stage_removed_mask = (raw_labels > 0) & ~(analysis_labels > 0)
        size_stage_display = labels_to_gray_display(
            analysis_labels,
            analysis_info,
        )
    else:
        raise ValueError(f"Unknown watershed mode: {watershed_mode}")

    return {
        "mask_cell": ~(analysis_labels > 0),
        "threshold_cell": threshold_value,
        "pre_watershed_cell_mask": pre_watershed_mask,
        "raw_labels": raw_labels,
        "raw_info": raw_info,
        "watershed_display": watershed_display,
        "size_filtered_labels": size_filtered_labels,
        "size_filtered_info": size_filtered_info,
        "size_removed_mask": size_removed_mask,
        "analysis_labels": analysis_labels,
        "analysis_info": analysis_info,
        "size_stage_display": size_stage_display,
        "stage_display": size_stage_display,
        "final_outline_reference_labels": raw_labels,
    }

def make_watershed_line_mask(labeled_mask_cell, line_width=5):
    if labeled_mask_cell is None or np.max(labeled_mask_cell) == 0:
        return np.zeros_like(labeled_mask_cell, dtype=bool)

    all_boundaries = find_boundaries(labeled_mask_cell, mode="inner")
    outer_boundaries = find_boundaries(labeled_mask_cell > 0, mode="inner")
    internal_boundaries = all_boundaries & ~outer_boundaries

    if line_width > 1:
        iterations = max(1, int(line_width) - 1)
        internal_boundaries = ndimage.binary_dilation(
            internal_boundaries,
            iterations=iterations
        )

        # Begrenzung auf Zellbereiche, damit keine Linien weit in den Hintergrund wachsen
        internal_boundaries = internal_boundaries & (labeled_mask_cell > 0)

    return internal_boundaries


def add_watershed_lines_to_plot(plot_mask, watershed_line_mask):
    plot_mask_with_lines = plot_mask.copy()

    if watershed_line_mask is None or not np.any(watershed_line_mask):
        return plot_mask_with_lines

    # In deinen Masken gilt: False = schwarze/aktuell sichtbare Zellfläche
    current_cell_pixels = ~plot_mask

    # Nur Linien behalten, die an aktuell sichtbare Zellflächen angrenzen
    line_support = ndimage.binary_dilation(current_cell_pixels, iterations=3)
    visible_line_mask = watershed_line_mask & line_support

    plot_mask_with_lines[visible_line_mask] = True  # True = weiße Linie im Graubild
    return plot_mask_with_lines


def threshold_method_label(threshold_setting):
    setting = str(threshold_setting).strip().lower()
    if setting == "mean":
        return "Mean"
    if setting == "otsu":
        return "Otsu"
    return "Manual"


def boolean_mask_to_display(mask):
    """
    Convert ViSCA-style boolean masks to display arrays.
    In most ViSCA masks, True = white/background and False = black/object.
    """
    return np.asarray(mask, dtype=float)


def outline_mask(mask, line_width=2):
    """Create a continuous object outline for display of excluded objects."""
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0 or not np.any(mask):
        return np.zeros_like(mask, dtype=bool)

    # Use filled and lightly closed objects so outlines are continuous rather than dotted.
    filled = ndimage.binary_fill_holes(mask)
    closed = ndimage.binary_closing(filled, iterations=2)

    # A small dilation draws the outline just outside the excluded object.
    # This keeps the excluded object interior white while the border remains closed.
    outer = ndimage.binary_dilation(closed, iterations=max(1, int(line_width)))
    inner = closed
    boundary = outer & ~inner
    return boundary


def labels_to_gray_display(label_image, watershed_segment_info=None, nucleus_holes_mask=None, excluded_outline_mask=None):
    """
    Display convention:
    - background: white
    - non-watershed objects: black
    - watershed-derived segments from the same original object: black/dark-to-light grey
    - excluded objects: white interior with continuous medium-grey outline
    """
    labels = np.asarray(label_image)
    display = np.ones(labels.shape, dtype=float)
    watershed_segment_info = watershed_segment_info or {}

    visible_labels = [int(x) for x in np.unique(labels) if x > 0]

    split_groups = {}
    for label in visible_labels:
        info = watershed_segment_info.get(label, {})
        raw_n = int(info.get("raw_number_of_segments_in_original_object", 1))
        if raw_n > 1:
            original_id = int(info.get("original_object_id", label))
            split_groups.setdefault(original_id, []).append(label)

    # The first watershed-derived segment may be black, but the first grey level
    # must be clearly distinct from black. The palette is capped below white.
    palette = [0.0, 0.45, 0.58, 0.70, 0.80, 0.88]
    label_to_shade = {}
    for original_id, group_labels in split_groups.items():
        for i, label in enumerate(sorted(group_labels)):
            label_to_shade[label] = palette[i % len(palette)]

    for label in visible_labels:
        shade = label_to_shade.get(label, 0.0)
        # Safety cap: object pixels are never displayed as pure white.
        display[labels == label] = min(float(shade), 0.88)

    if nucleus_holes_mask is not None:
        display[np.asarray(nucleus_holes_mask, dtype=bool)] = 1.0

    if excluded_outline_mask is not None:
        display[outline_mask(excluded_outline_mask, line_width=2)] = 0.30

    return display


def exclude_labels_from_watershed_split(labeled_split_labels, watershed_segment_info):
    """
    Build the Exclude-mode analysis mask from the same early watershed result
    used by Split mode. Any original object that was split by watershed
    (raw_number_of_segments_in_original_object > 1) is removed completely.
    Only not-watershed-affected objects continue through downstream filters.
    """
    out = np.zeros_like(labeled_split_labels, dtype=int)
    new_info = {}
    new_label = 1

    for region in measure.regionprops(labeled_split_labels):
        old_label = int(region.label)
        info = dict(watershed_segment_info.get(old_label, {}))
        raw_n = int(info.get("raw_number_of_segments_in_original_object", 1))

        if raw_n > 1:
            continue

        rows = region.coords[:, 0]
        cols = region.coords[:, 1]
        out[rows, cols] = new_label

        info["raw_number_of_segments_in_original_object"] = raw_n
        info["number_of_segments_in_original_object"] = 1
        new_info[new_label] = info
        new_label += 1

    return out, new_info


def filtered_labels_from_mask(label_image, keep_mask):
    out = np.asarray(label_image).copy()
    out[~np.asarray(keep_mask, dtype=bool)] = 0
    return out


def split_stage_from_cell_image(
    image_path,
    min_size_pixels,
    max_size_pixels,
    intensity_threshold_value,
    sigma,
    min_distance,
    h_value,
):
    """Run the Split-mode cellbody threshold + watershed stage."""
    return intensity_threshold_with_watershed_split(
        image_path,
        min_size_pixels,
        max_size_pixels,
        intensity_threshold_value,
        sigma=sigma,
        min_distance=min_distance,
        h_value=h_value,
    )


def get_label_neighbors(label_image):
    neighbors = {
        int(label): set()
        for label in np.unique(label_image)
        if label > 0
    }

    top = label_image[:-1, :]
    bottom = label_image[1:, :]

    vertical_mask = (
        (top != bottom)
        & (top > 0)
        & (bottom > 0)
    )

    vertical_pairs = np.stack(
        [top[vertical_mask], bottom[vertical_mask]],
        axis=1
    ) if np.any(vertical_mask) else np.empty((0, 2), dtype=int)

    left = label_image[:, :-1]
    right = label_image[:, 1:]

    horizontal_mask = (
        (left != right)
        & (left > 0)
        & (right > 0)
    )

    horizontal_pairs = np.stack(
        [left[horizontal_mask], right[horizontal_mask]],
        axis=1
    ) if np.any(horizontal_mask) else np.empty((0, 2), dtype=int)

    all_pairs = np.vstack([vertical_pairs, horizontal_pairs])

    if len(all_pairs) > 0:
        for label_a, label_b in np.unique(all_pairs, axis=0):
            label_a = int(label_a)
            label_b = int(label_b)

            neighbors.setdefault(label_a, set()).add(label_b)
            neighbors.setdefault(label_b, set()).add(label_a)

    return neighbors


def create_watershed_metadata_dataframe(
    base_image_name,
    labeled_mask_cell,
    labeled_unchanged_objects_mask,
    labeled_mask_nucleus,
    cell_dataframe,
    watershed_segment_info,
    scale,
):
    columns = [
        "Image Name",
        "Cell Number",
        "Cell Area (µm²)",
        "Nucleus Area (µm²)",
        "Watershed Segment ID",
        "Original Pre-Watershed Object ID",
        "Original Pre-Watershed Object Area (µm²)",
        "Segment Fraction Of Original Object",
        "Raw Number Of Watershed Segments In Original Object",
        "Original Object Was Split",
        "Number Of Segments In Original Object",
        "Segment Complete Nuclei Count",
        "Segment Has Complete Nucleus",
        "Number Of Nucleated Segments In Original Object",
        "Number Of Anucleated Segments In Original Object",
        "Border To Nucleated Segment",
        "Border To Anucleated Segment",
        "Split Group Type",
    ]

    if watershed_segment_info is None:
        return pd.DataFrame(columns=columns)

    if len(watershed_segment_info) == 0:
        return pd.DataFrame(columns=columns)

    nucleus_regions = {
        region.label: region
        for region in measure.regionprops(labeled_mask_nucleus)
    }

    # Count complete valid nuclei per watershed segment.
    segment_complete_nucleus_labels = {}

    for segment_region in measure.regionprops(labeled_mask_cell):
        segment_label = int(segment_region.label)
        segment_coords = segment_region.coords

        overlapping_nuclei = labeled_mask_nucleus[
            segment_coords[:, 0],
            segment_coords[:, 1],
        ]

        unique_nuclei = np.unique(overlapping_nuclei[overlapping_nuclei > 0])

        complete_nuclei = []

        for nucleus_label in unique_nuclei:
            nucleus_label = int(nucleus_label)
            nucleus_region = nucleus_regions.get(nucleus_label)

            if nucleus_region is None:
                continue

            nucleus_area_in_segment = np.sum(overlapping_nuclei == nucleus_label)

            if nucleus_area_in_segment == nucleus_region.area:
                complete_nuclei.append(nucleus_label)

        segment_complete_nucleus_labels[segment_label] = complete_nuclei

    original_to_segments = {}

    for segment_label, info in watershed_segment_info.items():
        original_id = info["original_object_id"]
        original_to_segments.setdefault(original_id, []).append(int(segment_label))

    neighbors = get_label_neighbors(labeled_mask_cell)

    rows = []

    for index, cell_region in enumerate(measure.regionprops(labeled_unchanged_objects_mask)):
        if index not in cell_dataframe.index:
            continue

        segment_label = int(cell_region.label)

        segment_info = watershed_segment_info.get(segment_label)

        if segment_info is None:
            continue

        original_id = segment_info["original_object_id"]
        segments_in_original = original_to_segments.get(original_id, [segment_label])

        n_segments_in_original = len(segments_in_original)

        original_object_area_pixels = segment_info.get("original_object_area_pixels", 0)
        original_object_area_um2 = original_object_area_pixels / (scale ** 2)

        cell_area_um2 = cell_dataframe.loc[index, "Cell Area (µm²)"]

        if original_object_area_um2 > 0:
            segment_fraction_of_original = cell_area_um2 / original_object_area_um2
        else:
            segment_fraction_of_original = np.nan

        n_nucleated_segments = sum(
            len(segment_complete_nucleus_labels.get(segment, [])) > 0
            for segment in segments_in_original
        )

        n_anucleated_segments = sum(
            len(segment_complete_nucleus_labels.get(segment, [])) == 0
            for segment in segments_in_original
        )

        neighboring_segments = [
            neighbor
            for neighbor in neighbors.get(segment_label, set())
            if watershed_segment_info.get(neighbor, {}).get("original_object_id") == original_id
        ]

        border_to_nucleated_segment = any(
            len(segment_complete_nucleus_labels.get(neighbor, [])) > 0
            for neighbor in neighboring_segments
        )

        border_to_anucleated_segment = any(
            len(segment_complete_nucleus_labels.get(neighbor, [])) == 0
            for neighbor in neighboring_segments
        )

        original_object_was_split = n_segments_in_original > 1

        if not original_object_was_split:
            split_group_type = "none"
        elif border_to_nucleated_segment and border_to_anucleated_segment:
            split_group_type = "mixed"
        elif border_to_nucleated_segment:
            split_group_type = "nucleated_neighbor"
        elif border_to_anucleated_segment:
            split_group_type = "anucleated_neighbor"
        else:
            split_group_type = "split_no_direct_remaining_neighbor"

        segment_complete_nuclei_count = len(
            segment_complete_nucleus_labels.get(segment_label, [])
        )

        rows.append(
            {
                "Image Name": base_image_name,
                "Cell Number": cell_dataframe.loc[index, "Cell Number"],
                "Cell Area (µm²)": cell_dataframe.loc[index, "Cell Area (µm²)"],
                "Nucleus Area (µm²)": cell_dataframe.loc[index, "Nucleus Area (µm²)"],
                "Watershed Segment ID": segment_label,
                "Original Pre-Watershed Object ID": original_id,
                "Original Pre-Watershed Object Area (µm²)": round(original_object_area_um2, 2),
                "Segment Fraction Of Original Object": round(segment_fraction_of_original, 4),
                "Raw Number Of Watershed Segments In Original Object": segment_info.get(
                    "raw_number_of_segments_in_original_object", 1
                ),
                "Original Object Was Split": "yes" if original_object_was_split else "no",
                "Number Of Segments In Original Object": n_segments_in_original,
                "Segment Complete Nuclei Count": segment_complete_nuclei_count,
                "Segment Has Complete Nucleus": "yes" if segment_complete_nuclei_count > 0 else "no",
                "Number Of Nucleated Segments In Original Object": n_nucleated_segments,
                "Number Of Anucleated Segments In Original Object": n_anucleated_segments,
                "Border To Nucleated Segment": "yes" if border_to_nucleated_segment else "no",
                "Border To Anucleated Segment": "yes" if border_to_anucleated_segment else "no",
                "Split Group Type": split_group_type,
            }
        )

    return pd.DataFrame(rows, columns=columns)

def mononuclear_filter(labeled_mask_cell, labeled_mask_nucleus):
    nucleus_overlap_mask = np.zeros_like(labeled_mask_cell, dtype=bool)
    nucleus_regions = {region.label: region for region in measure.regionprops(labeled_mask_nucleus)}

    for cell_region in measure.regionprops(labeled_mask_cell):
        cell_coords = cell_region.coords
        overlapping_nuclei = labeled_mask_nucleus[cell_coords[:, 0], cell_coords[:, 1]]
        unique_nuclei = np.unique(overlapping_nuclei[overlapping_nuclei > 0])

        if len(unique_nuclei) == 1:
            nucleus_label = unique_nuclei[0]
            nucleus_region = nucleus_regions[nucleus_label]
            nucleus_area_in_cell = np.sum(overlapping_nuclei == nucleus_label)

            if nucleus_area_in_cell == nucleus_region.area:
                nucleus_overlap_mask[cell_coords[:, 0], cell_coords[:, 1]] = True

    return nucleus_overlap_mask


def calculate_circularity(perimeter, area):
    if perimeter == 0:
        return 0
    return 4 * pi * area / (perimeter ** 2)


def circularity_filter(input_mask, labeled_mask_nucleus, labeled_mask_cell, circularity_threshold=0.3):
    circularity_mask = np.zeros_like(input_mask, dtype=bool)

    for region in measure.regionprops(labeled_mask_nucleus):
        circularity = calculate_circularity(region.perimeter, region.area)
        if 0 < circularity <= circularity_threshold:
            circularity_mask[region.coords[:, 0], region.coords[:, 1]] = True

    corrected_subtraction_mask = np.ones_like(input_mask, dtype=bool)
    corrected_subtraction_mask[input_mask == True] = False
    corrected_subtraction_mask[circularity_mask == True] = True

    circularity_overlap_mask = np.ones_like(corrected_subtraction_mask, dtype=bool)

    for region in measure.regionprops(labeled_mask_cell):
        mitotic_overlap = np.any(circularity_mask[region.coords[:, 0], region.coords[:, 1]])
        if not mitotic_overlap:
            circularity_overlap_mask[region.coords[:, 0], region.coords[:, 1]] = corrected_subtraction_mask[
                region.coords[:, 0], region.coords[:, 1]
            ]

    return circularity_mask, circularity_overlap_mask


def watershed_filter(input_mask, sigma=6, min_distance=30, h_value=3):
    distance = ndimage.distance_transform_edt(~input_mask)
    smoothed_distance = gaussian(distance, sigma=sigma)
    h_maxima_transformed = h_maxima(smoothed_distance, h=h_value)

    local_maxi = peak_local_max(
        h_maxima_transformed,
        footprint=np.ones((3, 3)),
        min_distance=min_distance,
        exclude_border=False,
    )

    markers = np.zeros_like(input_mask, dtype=int)
    if len(local_maxi) > 0:
        markers[tuple(local_maxi.T)] = np.arange(1, len(local_maxi) + 1)

    watershed_labels = watershed(-smoothed_distance, markers, mask=~input_mask)

    unchanged_objects_mask = np.zeros_like(watershed_labels, dtype=bool)
    pre_watershed_labels = measure.label(~input_mask)

    for region in measure.regionprops(pre_watershed_labels):
        watershed_labels_in_region = watershed_labels[region.coords[:, 0], region.coords[:, 1]]
        unique_labels = np.unique(watershed_labels_in_region[watershed_labels_in_region > 0])

        if len(unique_labels) == 1:
            for coord in region.coords:
                unchanged_objects_mask[coord[0], coord[1]] = True

    return unchanged_objects_mask


def measure_cells(input_mask, labeled_mask_nucleus, scale, ratio_threshold):
    cell_sizes = []
    nucleus_sizes = []

    nucleus_regions = {region.label: region for region in measure.regionprops(labeled_mask_nucleus)}

    for cell_region in measure.regionprops(input_mask):
        cell_size = cell_region.area / (scale ** 2)
        cell_coords = cell_region.coords
        overlapping_nuclei = labeled_mask_nucleus[cell_coords[:, 0], cell_coords[:, 1]]
        unique_nuclei = np.unique(overlapping_nuclei[overlapping_nuclei > 0])

        if len(unique_nuclei) == 1:
            nucleus_label = unique_nuclei[0]
            nucleus_region = nucleus_regions[nucleus_label]
            nucleus_size = nucleus_region.area / (scale ** 2)

            cell_sizes.append(cell_size)
            nucleus_sizes.append(nucleus_size)

    cell_dataframe = pd.DataFrame(
        {
            "Cell Number": range(1, len(cell_sizes) + 1),
            "Cell Area (µm²)": cell_sizes,
            "Nucleus Area (µm²)": nucleus_sizes,
        }
    )

    cell_dataframe["Cell Area (µm²)"] = cell_dataframe["Cell Area (µm²)"].round(2)
    cell_dataframe["Nucleus Area (µm²)"] = cell_dataframe["Nucleus Area (µm²)"].round(2)

    cell_dataframe = cell_dataframe[
        cell_dataframe["Nucleus Area (µm²)"] / cell_dataframe["Cell Area (µm²)"] < ratio_threshold
    ]

    cell_dataframe["Cell Number"] = range(1, len(cell_dataframe) + 1)
    return cell_dataframe


def circularity_filter_preview(input_mask, circularity_threshold=0.3):
    input_mask = ~input_mask
    labeled_mask = measure.label(input_mask)
    final_mask = np.zeros_like(input_mask, dtype=bool)

    for region in measure.regionprops(labeled_mask):
        circularity = calculate_circularity(region.perimeter, region.area)
        if circularity > circularity_threshold:
            final_mask[region.coords[:, 0], region.coords[:, 1]] = True

    return ~final_mask


def ratio_filter_preview(
    cell_image_path,
    nucleus_image_path,
    scale,
    intensity_threshold_cell_value,
    min_size_cell,
    max_size_cell,
    intensity_threshold_nucleus_value,
    min_size_nucleus,
    max_size_nucleus,
    circularity_threshold,
    watershed_mode,
    watershed_sigma,
    watershed_min_distance,
    watershed_h_value,
):
    min_cell_px = min_size_cell * (scale ** 2)
    max_cell_px = max_size_cell * (scale ** 2)
    min_nucleus_px = min_size_nucleus * (scale ** 2)
    max_nucleus_px = max_size_nucleus * (scale ** 2)

    cell_stage = run_cellbody_stage(
        cell_image_path,
        min_cell_px,
        max_cell_px,
        intensity_threshold_cell_value,
        watershed_mode=watershed_mode,
        sigma=watershed_sigma,
        min_distance=watershed_min_distance,
        h_value=watershed_h_value,
    )

    labeled_mask_cell = cell_stage["analysis_labels"]
    watershed_segment_info = cell_stage["analysis_info"]

    mask_nucleus, _ = intensity_threshold(
        nucleus_image_path,
        min_nucleus_px,
        max_nucleus_px,
        intensity_threshold_nucleus_value,
    )

    labeled_mask_nucleus = measure.label(~mask_nucleus)

    nucleus_overlap_mask = mononuclear_filter(labeled_mask_cell, labeled_mask_nucleus)

    circularity_mask, circularity_overlap_mask = circularity_filter(
        nucleus_overlap_mask,
        labeled_mask_nucleus,
        labeled_mask_cell,
        circularity_threshold=circularity_threshold,
    )

    valid_cell_mask = ~circularity_overlap_mask
    if watershed_mode in {"Split", "Exclude"}:
        labeled_valid_cells = labeled_mask_cell.copy()
        labeled_valid_cells[~valid_cell_mask] = 0
    else:
        labeled_valid_cells = measure.label(valid_cell_mask)

    nucleus_regions = {
        region.label: region
        for region in measure.regionprops(labeled_mask_nucleus)
    }

    valid_objects = []

    for region in measure.regionprops(labeled_valid_cells):
        cell_coords = region.coords
        overlapping_nuclei = labeled_mask_nucleus[cell_coords[:, 0], cell_coords[:, 1]]
        unique_nuclei = np.unique(overlapping_nuclei[overlapping_nuclei > 0])

        if len(unique_nuclei) != 1:
            continue

        nucleus_label = unique_nuclei[0]
        nucleus_region = nucleus_regions.get(nucleus_label)

        if nucleus_region is None:
            continue

        nucleus_area_in_cell = np.sum(overlapping_nuclei == nucleus_label)
        if nucleus_area_in_cell != nucleus_region.area:
            continue

        cell_area = region.area / (scale ** 2)
        nucleus_area = nucleus_region.area / (scale ** 2)

        valid_objects.append(
            {
                "cell_coords": cell_coords,
                "nucleus_coords": nucleus_region.coords,
                "cell_area": cell_area,
                "nucleus_area": nucleus_area,
                "label": int(region.label),
            }
        )

    preview_shape = cell_stage["pre_watershed_cell_mask"].shape
    return valid_objects, preview_shape, labeled_valid_cells, watershed_segment_info



# Output functions

def save_csv(cell_dataframe, target_dir, base_image_name, watershed_mode=None):
    cell_dataframe = cell_dataframe.copy()
    cell_dataframe.insert(0, "Image Name", base_image_name)
    if watershed_mode is not None and "Watershed Mode" not in cell_dataframe.columns:
        cell_dataframe.insert(1, "Watershed Mode", watershed_mode)

    csv_subfolder = os.path.join(target_dir, "Area_csv")
    ensure_directory(csv_subfolder)

    csv_filename = f"{base_image_name}_cell_areas.csv"
    csv_file_path = os.path.join(csv_subfolder, csv_filename)

    cell_dataframe.to_csv(csv_file_path, index=False, sep=";", decimal=",")
    print(f"Saved table to {csv_file_path}")




def create_non_circular_nuclei_dataframe(
    base_image_name,
    watershed_mode,
    mononuclear_labels,
    labeled_mask_nucleus,
    circularity_mask,
    scale,
    circularity_threshold,
):
    """Create per-cell output for cells excluded by the nucleus circularity filter."""
    columns = [
        "Image Name",
        "Watershed Mode",
        "Non-circular Cell ID",
        "Cell Area (µm²)",
        "Nucleus Area (µm²)",
        "Nucleus Circularity",
        "Circularity Threshold",
    ]

    if mononuclear_labels is None or np.max(mononuclear_labels) == 0:
        return pd.DataFrame(columns=columns)

    nucleus_regions = {
        int(region.label): region
        for region in measure.regionprops(labeled_mask_nucleus)
    }

    rows = []
    counter = 1

    for cell_region in measure.regionprops(mononuclear_labels):
        cell_coords = cell_region.coords
        rr = cell_coords[:, 0]
        cc = cell_coords[:, 1]

        # Only cells that are removed at this step because they overlap a non-circular nucleus.
        if not np.any(circularity_mask[rr, cc]):
            continue

        overlapping_nuclei = labeled_mask_nucleus[rr, cc]
        unique_nuclei = np.unique(overlapping_nuclei[overlapping_nuclei > 0])

        complete_nuclei = []
        for nucleus_label in unique_nuclei:
            nucleus_label = int(nucleus_label)
            nucleus_region = nucleus_regions.get(nucleus_label)
            if nucleus_region is None:
                continue

            nucleus_area_in_cell = np.sum(overlapping_nuclei == nucleus_label)
            if nucleus_area_in_cell == nucleus_region.area:
                circularity = calculate_circularity(nucleus_region.perimeter, nucleus_region.area)
                complete_nuclei.append((nucleus_label, nucleus_region, circularity))

        non_circular_complete_nuclei = [
            item for item in complete_nuclei
            if 0 < item[2] <= circularity_threshold
        ]

        if non_circular_complete_nuclei:
            nucleus_label, nucleus_region, nucleus_circularity = sorted(
                non_circular_complete_nuclei, key=lambda x: x[2]
            )[0]
        elif complete_nuclei:
            # Fallback: should rarely be needed, but keeps the output robust if masks disagree.
            nucleus_label, nucleus_region, nucleus_circularity = sorted(
                complete_nuclei, key=lambda x: x[2]
            )[0]
        else:
            continue

        rows.append(
            {
                "Image Name": base_image_name,
                "Watershed Mode": watershed_mode,
                "Non-circular Cell ID": f"non-circular#{counter}",
                "Cell Area (µm²)": round(cell_region.area / (scale ** 2), 2),
                "Nucleus Area (µm²)": round(nucleus_region.area / (scale ** 2), 2),
                "Nucleus Circularity": round(float(nucleus_circularity), 4),
                "Circularity Threshold": circularity_threshold,
            }
        )
        counter += 1

    return pd.DataFrame(rows, columns=columns)


def save_non_circular_nuclei_csv(non_circular_dataframe, target_dir, base_image_name):
    csv_subfolder = os.path.join(target_dir, "Non_circular_nuclei_csv")
    ensure_directory(csv_subfolder)

    csv_filename = f"{base_image_name}_non_circular_nuclei.csv"
    csv_file_path = os.path.join(csv_subfolder, csv_filename)

    non_circular_dataframe.to_csv(csv_file_path, index=False, sep=";", decimal=",")
    print(f"Saved non-circular nuclei table to {csv_file_path}")

def save_watershed_metadata_csv(metadata_dataframe, target_dir, base_image_name):
    metadata_subfolder = os.path.join(target_dir, "Watershed_metadata_csv")
    ensure_directory(metadata_subfolder)

    metadata_filename = f"{base_image_name}_watershed_metadata.csv"
    metadata_file_path = os.path.join(metadata_subfolder, metadata_filename)

    metadata_dataframe.to_csv(
        metadata_file_path,
        index=False,
        sep=";",
        decimal=",",
    )

    print(f"Saved watershed metadata to {metadata_file_path}")


def save_plots(plots, cell_dataframe, final_label_image, target_dir, base_image_name):
    png_subfolder = os.path.join(target_dir, "Masks_png")
    ensure_directory(png_subfolder)

    png_filename = f"{base_image_name}_mask.png"
    png_file_path = os.path.join(png_subfolder, png_filename)

    plt.figure(figsize=(17, 7.8))

    plot_cols = 4
    plot_rows = 2

    items = list(plots.items())
    for index, (title, plot) in enumerate(items[:8]):
        plt.subplot(plot_rows, plot_cols, index + 1)
        plt.imshow(plot, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        plt.title(title)
        plt.axis("on")

        if index == 7 and final_label_image is not None:
            for source_index, region in enumerate(measure.regionprops(final_label_image)):
                if source_index in cell_dataframe.index:
                    cell_num = int(cell_dataframe.loc[source_index, "Cell Number"])
                    minr, minc, maxr, maxc = region.bbox
                    label_x = maxc + 5
                    label_y = (minr + maxr) / 2
                    plt.text(
                        label_x,
                        label_y,
                        f"{cell_num}",
                        color=required_color,
                        fontsize=8,
                        ha="left",
                        va="center",
                    )

    plt.tight_layout(w_pad=0.4, h_pad=1.0)
    plt.savefig(png_file_path, dpi=200)
    plt.close()

    print(f"Saved figure to {png_file_path}")


def save_parameters_to_folders(cell_folders, parameters):
    for cell_folder in cell_folders:
        parameter_file = os.path.join(cell_folder, "parameters.txt")
        with open(parameter_file, "w", encoding="utf-8") as f:
            for key, value in parameters.items():
                f.write(f"{key}: {value}\n")


def aggregate_csvs(cell_folders):
    for cell_folder in cell_folders:
        csv_path = os.path.join(cell_folder, "Area_csv")
        if not os.path.isdir(csv_path):
            continue

        all_files = os.listdir(csv_path)
        all_csv_files = sorted([file for file in all_files if file.endswith(".csv") and file != "aggregated_data.csv"])

        if not all_csv_files:
            continue

        all_csv_dfs = [
            pd.read_csv(os.path.join(csv_path, file), sep=";", decimal=",")
            for file in all_csv_files
        ]
        all_csv_df = pd.concat(all_csv_dfs, ignore_index=True)
        all_csv_df_sorted = all_csv_df.sort_values(by=["Image Name", "Cell Number"])

        all_csv_df_sorted.to_csv(
            os.path.join(csv_path, "aggregated_data.csv"),
            index=False,
            sep=";",
            decimal=",",
        )



def aggregate_non_circular_nuclei_csvs(cell_folders):
    for cell_folder in cell_folders:
        csv_path = os.path.join(cell_folder, "Non_circular_nuclei_csv")
        if not os.path.isdir(csv_path):
            continue

        all_files = os.listdir(csv_path)
        all_csv_files = sorted([
            file for file in all_files
            if file.endswith(".csv") and file != "aggregated_non_circular_nuclei.csv"
        ])

        if not all_csv_files:
            continue

        all_csv_dfs = [
            pd.read_csv(os.path.join(csv_path, file), sep=";", decimal=",")
            for file in all_csv_files
        ]
        all_csv_df = pd.concat(all_csv_dfs, ignore_index=True)

        if len(all_csv_df) > 0 and {"Image Name", "Non-circular Cell ID"}.issubset(all_csv_df.columns):
            all_csv_df = all_csv_df.sort_values(by=["Image Name", "Non-circular Cell ID"])

        all_csv_df.to_csv(
            os.path.join(csv_path, "aggregated_non_circular_nuclei.csv"),
            index=False,
            sep=";",
            decimal=",",
        )

        print(
            f"Saved aggregated non-circular nuclei table to "
            f"{os.path.join(csv_path, 'aggregated_non_circular_nuclei.csv')}"
        )

def aggregate_watershed_metadata_csvs(cell_folders):
    for cell_folder in cell_folders:
        metadata_path = os.path.join(cell_folder, "Watershed_metadata_csv")

        if not os.path.isdir(metadata_path):
            continue

        all_files = os.listdir(metadata_path)

        all_metadata_csv_files = sorted(
            [
                file
                for file in all_files
                if file.endswith(".csv")
                and file != "aggregated_watershed_metadata.csv"
            ]
        )

        if not all_metadata_csv_files:
            continue

        all_metadata_dfs = [
            pd.read_csv(
                os.path.join(metadata_path, file),
                sep=";",
                decimal=",",
            )
            for file in all_metadata_csv_files
        ]

        all_metadata_df = pd.concat(all_metadata_dfs, ignore_index=True)

        all_metadata_df_sorted = all_metadata_df.sort_values(
            by=["Image Name", "Cell Number"]
        )

        all_metadata_df_sorted.to_csv(
            os.path.join(metadata_path, "aggregated_watershed_metadata.csv"),
            index=False,
            sep=";",
            decimal=",",
        )

        print(
            f"Saved aggregated watershed metadata to "
            f"{os.path.join(metadata_path, 'aggregated_watershed_metadata.csv')}"
        )

# Processing pipeline

def process_image_pair(image_path_cell, image_path_nucleus, parameters):
    scale = parameters["scale"]
    intensity_threshold_cell_value = parameters["intensity_threshold_cell"]
    min_size_cell = parameters["min_size_cell"]
    max_size_cell = parameters["max_size_cell"]
    intensity_threshold_nucleus_value = parameters["intensity_threshold_nucleus"]
    min_size_nucleus = parameters["min_size_nucleus"]
    max_size_nucleus = parameters["max_size_nucleus"]
    circularity_threshold = parameters["circularity_threshold"]
    ratio_threshold = parameters["ratio_threshold"]
    watershed_mode = parameters.get("watershed_mode", "No Watershed")
    watershed_sigma = parameters["watershed_sigma"]
    watershed_min_distance = parameters["watershed_min_distance"]
    watershed_h_value = parameters["watershed_h_value"]

    source_dir = os.path.dirname(image_path_cell)
    min_cell_px = min_size_cell * (scale ** 2)
    max_cell_px = max_size_cell * (scale ** 2)
    min_nucleus_px = min_size_nucleus * (scale ** 2)
    max_nucleus_px = max_size_nucleus * (scale ** 2)

    cell_stage = run_cellbody_stage(
        image_path_cell,
        min_cell_px,
        max_cell_px,
        intensity_threshold_cell_value,
        watershed_mode=watershed_mode,
        sigma=watershed_sigma,
        min_distance=watershed_min_distance,
        h_value=watershed_h_value,
    )

    threshold_cell = cell_stage["threshold_cell"]
    pre_watershed_cell_mask = cell_stage["pre_watershed_cell_mask"]
    labeled_mask_cell = cell_stage["analysis_labels"]
    watershed_segment_info = cell_stage["analysis_info"]
    watershed_display = cell_stage["watershed_display"]
    size_stage_display = cell_stage["size_stage_display"]
    final_outline_reference_labels = cell_stage["final_outline_reference_labels"]

    if watershed_mode == "No Watershed":
        watershed_stage_title = "No Watershed applied"
    elif watershed_mode == "Split":
        watershed_stage_title = (
            "Watershed split mode "
            f"(σ={watershed_title_number(watershed_sigma)}, "
            f"d={watershed_title_number(watershed_min_distance)}, "
            f"h={watershed_title_number(watershed_h_value)})"
        )
    else:
        watershed_stage_title = (
            "Watershed exclude mode "
            f"(σ={watershed_title_number(watershed_sigma)}, "
            f"d={watershed_title_number(watershed_min_distance)}, "
            f"h={watershed_title_number(watershed_h_value)})"
        )

    size_stage_title = cell_size_range_title(min_size_cell, max_size_cell)
    if watershed_mode == "Exclude":
        size_stage_title = (
            f"{size_stage_title}\n"
            "& Watershed-affected objects excluded"
        )

    mask_nucleus, threshold_nucleus = intensity_threshold(
        image_path_nucleus,
        min_nucleus_px,
        max_nucleus_px,
        intensity_threshold_nucleus_value,
    )

    labeled_mask_nucleus = measure.label(~mask_nucleus)

    nucleus_overlap_mask = mononuclear_filter(labeled_mask_cell, labeled_mask_nucleus)
    mononuclear_labels = filtered_labels_from_mask(labeled_mask_cell, nucleus_overlap_mask)
    mononuclear_display = labels_to_gray_display(mononuclear_labels, watershed_segment_info)

    circularity_mask, circularity_overlap_mask = circularity_filter(
        nucleus_overlap_mask,
        labeled_mask_nucleus,
        labeled_mask_cell,
        circularity_threshold=circularity_threshold,
    )

    circular_valid_mask = ~circularity_overlap_mask
    circular_labels = filtered_labels_from_mask(labeled_mask_cell, circular_valid_mask)

    # Display for the circularity step:
    # valid cells remain black/grey; cells excluded because of a non-circular nucleus
    # are shown as white objects with grey outline and the non-circular nucleus in black.
    circular_rejected_mask = np.zeros_like(circular_valid_mask, dtype=bool)
    circular_display_valid_labels = np.zeros_like(labeled_mask_cell, dtype=int)
    for region in measure.regionprops(mononuclear_labels):
        rows = region.coords[:, 0]
        cols = region.coords[:, 1]
        if np.any(circularity_mask[rows, cols]):
            circular_rejected_mask[rows, cols] = True
        else:
            circular_display_valid_labels[rows, cols] = int(region.label)

    circular_cells_display = labels_to_gray_display(
        circular_display_valid_labels,
        watershed_segment_info,
        excluded_outline_mask=circular_rejected_mask,
    )
    circular_cells_display[circularity_mask & circular_rejected_mask] = 0.0

    unchanged_objects_mask = circular_valid_mask
    if watershed_mode in {"Split", "Exclude"}:
        labeled_unchanged_objects_mask = labeled_mask_cell.copy()
        labeled_unchanged_objects_mask[~unchanged_objects_mask] = 0
    else:
        labeled_unchanged_objects_mask = measure.label(unchanged_objects_mask)

    cell_dataframe = measure_cells(
        labeled_unchanged_objects_mask,
        labeled_mask_nucleus,
        scale,
        ratio_threshold,
    )

    final_mask = np.zeros_like(labeled_unchanged_objects_mask, dtype=bool)
    final_nucleus_holes_mask = np.zeros_like(labeled_unchanged_objects_mask, dtype=bool)

    for index, region in enumerate(measure.regionprops(labeled_unchanged_objects_mask)):
        if index in cell_dataframe.index:
            final_mask[region.coords[:, 0], region.coords[:, 1]] = True

            cell_coords = region.coords
            overlapping_nuclei = labeled_mask_nucleus[cell_coords[:, 0], cell_coords[:, 1]]
            unique_nuclei = np.unique(overlapping_nuclei[overlapping_nuclei > 0])

            if len(unique_nuclei) == 1:
                nucleus_label = unique_nuclei[0]
                nucleus_region_mask = labeled_mask_nucleus == nucleus_label
                final_nucleus_holes_mask[nucleus_region_mask] = True

    ratio_labels = labeled_unchanged_objects_mask.copy()
    ratio_labels[~final_mask] = 0
    excluded_mask = (final_outline_reference_labels > 0) & ~final_mask
    ratio_display = labels_to_gray_display(
        ratio_labels,
        watershed_segment_info,
        nucleus_holes_mask=final_nucleus_holes_mask,
    )
    final_display = labels_to_gray_display(
        ratio_labels,
        watershed_segment_info,
        nucleus_holes_mask=final_nucleus_holes_mask,
        excluded_outline_mask=excluded_mask,
    )

    final_label_image = ratio_labels

    base_image_name = os.path.splitext(os.path.basename(image_path_cell))[0]

    non_circular_dataframe = create_non_circular_nuclei_dataframe(
        base_image_name,
        watershed_mode,
        mononuclear_labels,
        labeled_mask_nucleus,
        circularity_mask,
        scale,
        circularity_threshold,
    )
    save_non_circular_nuclei_csv(
        non_circular_dataframe,
        source_dir,
        base_image_name,
    )

    save_csv(cell_dataframe, source_dir, base_image_name, watershed_mode=watershed_mode)

    cell_method = threshold_method_label(intensity_threshold_cell_value)
    nucleus_method = threshold_method_label(intensity_threshold_nucleus_value)

    mononuclear_title = "Mononuclear cells"

    plots = {
        f"Cell intensity threshold\n({cell_method}; > {threshold_cell:.0f})": boolean_mask_to_display(pre_watershed_cell_mask),
        watershed_stage_title: watershed_display,
        size_stage_title: size_stage_display,
        f"Nucleus intensity threshold\n({nucleus_method}; > {threshold_nucleus:.0f})": boolean_mask_to_display(mask_nucleus),
        mononuclear_title: mononuclear_display,
        f"Nucleus circularity ≤ {circularity_threshold} excluded": circular_cells_display,
        f"Nucleus / Cellbody Area Ratio < {ratio_threshold}": ratio_display,
        "Final mask": final_display,
    }

    save_plots(
        plots,
        cell_dataframe,
        labeled_unchanged_objects_mask,
        source_dir,
        base_image_name,
    )

def process_single_image_pair(args):
    load_analysis_modules()

    image_path_cell, image_path_nucleus, parameters = args
    process_image_pair(image_path_cell, image_path_nucleus, parameters)

    return 1

def run_batch_processing(image_paths_pairs, parameters):
    if not image_paths_pairs:
        raise ValueError("No image pairs were found.")

    progress_root = tk.Tk()
    progress_root.title("Processing Images")
    progress_root.resizable(False, False)

    total_images = len(image_paths_pairs)
    completed_images = 0
    pool_closed = False

    progress = ttk.Progressbar(
        progress_root,
        length=340,
        mode="determinate",
        maximum=total_images,
        value=0
    )
    progress.pack(pady=(12, 6), padx=14)

    status_label = ttk.Label(
        progress_root,
        text=f"Processing images... 0 / {total_images} (0%)"
    )
    status_label.pack(pady=(0, 12))

    def update_progress():
        percent = int((completed_images / total_images) * 100)
        progress["value"] = completed_images
        status_label.config(
            text=f"Processing images... {completed_images} / {total_images} ({percent}%)"
        )
        progress_root.update_idletasks()

    def on_image_done(_result):
        nonlocal completed_images, pool_closed

        completed_images += 1
        update_progress()

        if completed_images >= total_images and not pool_closed:
            pool_closed = True
            pool.close()

            def finish():
                pool.join()
                progress["value"] = total_images
                status_label.config(
                    text=f"Processing images... {total_images} / {total_images} (100%)"
                )
                progress_root.update_idletasks()
                progress_root.after(500, progress_root.destroy)
                print("All images processed")

            progress_root.after(100, finish)

    def on_error(error):
        nonlocal pool_closed

        if not pool_closed:
            pool_closed = True
            pool.terminate()

        progress_root.destroy()
        messagebox.showerror("Processing Error", str(error))

    tasks = [
        (image_path_cell, image_path_nucleus, parameters)
        for image_path_cell, image_path_nucleus in image_paths_pairs
    ]

    pool = multiprocessing.Pool(processes=NUM_PROCESSES)

    for task in tasks:
        pool.apply_async(
            process_single_image_pair,
            args=(task,),
            callback=on_image_done,
            error_callback=on_error
        )

    update_progress()
    progress_root.mainloop()


# Folder selection

class DirectorySelectorTK:
    def __init__(self, title, folder_label):
        self.selected_directories = []
        self.root = tk.Tk()
        self.root.title(title)
        self.root.geometry("700x420")
        self.root.minsize(700, 420)

        main_frame = ttk.Frame(self.root, padding="12")
        main_frame.pack(fill="both", expand=True)

        heading = ttk.Label(main_frame, text=title, font=("Arial", 14, "bold"))
        heading.pack(anchor="w", pady=(0, 10))

        info_text = (
            f"Add one or more folders containing the {folder_label} images.\n"
            "Use 'Add folder' repeatedly if your images are distributed across multiple folders."
        )
        info_label = ttk.Label(main_frame, text=info_text, justify="left")
        info_label.pack(anchor="w", pady=(0, 10))

        self.listbox = tk.Listbox(main_frame, selectmode=tk.EXTENDED, width=100, height=14)
        self.listbox.pack(fill="both", expand=True, pady=(0, 10))

        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x")

        ttk.Button(button_frame, text="Add folder", command=self.add_folder).pack(side="left", padx=(0, 8))
        ttk.Button(button_frame, text="Remove selected", command=self.remove_selected).pack(side="left", padx=(0, 8))
        ttk.Button(button_frame, text="Done", command=self.finish).pack(side="right")

        self.root.protocol("WM_DELETE_WINDOW", self.close_without_selection)

    def add_folder(self):
        folder = filedialog.askdirectory(title="Select folder")
        if folder and folder not in self.selected_directories:
            self.selected_directories.append(folder)
            self.listbox.insert(tk.END, folder)

    def remove_selected(self):
        selected_indices = list(self.listbox.curselection())
        selected_indices.reverse()

        for index in selected_indices:
            self.listbox.delete(index)
            del self.selected_directories[index]

    def finish(self):
        if not self.selected_directories:
            messagebox.showwarning("No Selection", "No directories were selected.")
            return
        self.root.destroy()

    def close_without_selection(self):
        self.selected_directories = []
        self.root.destroy()

    def run(self):
        self.root.mainloop()
        return sorted(self.selected_directories)


def select_directories_tk(itype):
    selector = DirectorySelectorTK(
        title=f"Select the folder(s) containing the {itype} images",
        folder_label=itype,
    )
    return selector.run()


# Preview windows

class BasePreviewWindow:
    WINDOW_WIDTH = 430
    IMAGE_SIZE = 400
    SLIDER_LENGTH = 400
    MARGIN = 8
    SLIDER_PAD_Y = 4

    def setup_window(self, title, height):
        self.master.title(title)

        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        height = min(int(height), max(520, int(screen_height) - 80))
        x = screen_width // 2 - self.WINDOW_WIDTH // 2
        y = max(20, screen_height // 2 - height // 2)
        self.master.geometry(f"{self.WINDOW_WIDTH}x{height}+{x}+{y}")

    def create_image_label(self):
        image_label = tk.Label(self.master)
        image_label.pack(pady=self.MARGIN)
        return image_label

    def create_blank_preview(self):
        image = Image.new("L", (self.IMAGE_SIZE, self.IMAGE_SIZE), 255)
        image.thumbnail((self.IMAGE_SIZE, self.IMAGE_SIZE))
        photo = ImageTk.PhotoImage(image)
        self.image_label.config(image=photo)
        self.image_label.image = photo

    def set_preview_image(self, image_array):
        arr = np.asarray(image_array)
        if arr.dtype == bool:
            arr_uint8 = arr.astype(np.uint8) * 255
        else:
            arr_uint8 = (np.clip(arr.astype(float), 0, 1) * 255).astype(np.uint8)
        image = Image.fromarray(arr_uint8)
        image.thumbnail((self.IMAGE_SIZE, self.IMAGE_SIZE))
        photo = ImageTk.PhotoImage(image)
        self.image_label.config(image=photo)
        self.image_label.image = photo

    def create_slider(self, **kwargs):
        slider = tk.Scale(
            self.master,
            length=self.SLIDER_LENGTH,
            orient="horizontal",
            **kwargs,
        )
        slider.pack(pady=self.SLIDER_PAD_Y)
        enable_arrow_key_step(slider)
        return slider


class ParameterPreviewWindow(BasePreviewWindow):
    def __init__(
        self,
        master,
        scale,
        default_threshold,
        default_min,
        default_max,
        auto_threshold_method="otsu",
        initial_file_path=None,
        watershed_mode=None,
        watershed_sigma=6.0,
        watershed_min_distance=30,
        watershed_h_value=3.0,
        show_size_sliders=True,
    ):
        self.image = None
        self.master = master
        self.show_size_sliders = bool(show_size_sliders)
        self.setup_window(
            "Image Viewer",
            height=690 if self.show_size_sliders else 640,
        )

        self.scale = scale
        self.intensity_threshold = default_threshold
        self.min_size_threshold = default_min
        self.max_size_threshold = default_max
        self.auto_threshold_method = auto_threshold_method
        self.file_path = initial_file_path
        self.watershed_mode = watershed_mode
        self.watershed_sigma = watershed_sigma
        self.watershed_min_distance = watershed_min_distance
        self.watershed_h_value = watershed_h_value

        self.image_label = self.create_image_label()

        self.load_button = ttk.Button(self.master, text="Load Image", command=self.load_image)
        self.load_button.pack(pady=5)

        self.image_min = 0
        self.image_max = 1
        self.resolution = 1

        self.checked = tk.IntVar()
        self.checkbox = tk.Checkbutton(
            self.master,
            text="manual adjustment of intensity threshold",
            variable=self.checked,
            command=self.apply_checkbox,
        )
        self.checkbox.pack(pady=5)

        self.intensity_threshold_slider = self.create_slider(
            from_=self.image_min,
            to=self.image_max,
            resolution=self.resolution,
            command=self.update_parameter,
        )

        self.min_size_threshold_slider = None
        self.max_size_threshold_slider = None
        if self.show_size_sliders:
            self.min_size_threshold_slider = self.create_slider(
                from_=0,
                to=54,
                resolution=1,
                command=self.update_parameter,
            )
            self.min_size_threshold_slider.set(slider_position_value(self.min_size_threshold))

            self.max_size_threshold_slider = self.create_slider(
                from_=0,
                to=54,
                resolution=1,
                command=self.update_parameter,
            )
            self.max_size_threshold_slider.set(slider_position_value(self.max_size_threshold))

            self.min_size_threshold_slider.config(showvalue=0)
            self.max_size_threshold_slider.config(showvalue=0)

        self.intensity_threshold_slider.config(showvalue=0)

        self.save_button = ttk.Button(self.master, text="Save Threshold", command=self.save_threshold)
        self.save_button.pack(pady=5)

        self.create_blank_preview()

        try:
            int(default_threshold)
            self.checkbox.select()
        except ValueError:
            pass

        self.checkbox.config(state="disabled")
        self.intensity_threshold_slider.config(state="disabled")
        self.save_button.config(state="disabled")

        if self.show_size_sliders:
            self.min_size_threshold_slider.config(state="disabled")
            self.max_size_threshold_slider.config(state="disabled")

        self.fg = self.intensity_threshold_slider.cget("fg")
        self.troughcolor = self.intensity_threshold_slider.cget("troughcolor")

        self.intensity_threshold_slider.config(troughcolor="#f0f0f0")
        self.checkbox.config(fg="#808080")
        self.intensity_threshold_slider.config(fg="#808080")

        if self.show_size_sliders:
            self.min_size_threshold_slider.config(troughcolor="#f0f0f0")
            self.max_size_threshold_slider.config(troughcolor="#f0f0f0")
            self.min_size_threshold_slider.config(fg="#808080")
            self.max_size_threshold_slider.config(fg="#808080")

        if self.file_path:
            self.load_image_from_path(self.file_path)

    def save_threshold(self):
        self.master.destroy()

    def apply_checkbox(self):
        if self.checked.get():
            self.intensity_threshold_slider.config(state="normal")
            self.intensity_threshold_slider.config(troughcolor=self.troughcolor)
            self.intensity_threshold_slider.config(fg=self.fg)
            self.intensity_threshold = self.intensity_threshold_slider.get()
            self.update_parameter(0)
        else:
            self.intensity_threshold_slider.config(state="disabled")
            self.intensity_threshold_slider.config(troughcolor="#f0f0f0")
            self.intensity_threshold_slider.config(fg="#808080")
            self.intensity_threshold = self.auto_threshold_method
            if self.file_path:
                self.show_image(self.file_path)
            self.update_parameter(0)

    def update_parameter(self, _value):
        if self.intensity_threshold_slider.cget("state") != "disabled":
            self.intensity_threshold = self.intensity_threshold_slider.get()

        self.intensity_threshold_slider.config(label=f"Intensity Threshold: {self.intensity_threshold}")

        if self.show_size_sliders:
            self.min_size_threshold = stepped_slider_sequence(self.min_size_threshold_slider.get())
            self.max_size_threshold = stepped_slider_sequence(self.max_size_threshold_slider.get())
            self.min_size_threshold_slider.config(label=f"Min Size Threshold: {self.min_size_threshold:.0f} µm²")
            self.max_size_threshold_slider.config(label=f"Max Size Threshold: {self.max_size_threshold:.0f} µm²")

        if self.file_path:
            self.show_image(self.file_path)

    def load_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.load_image_from_path(file_path)

    def load_image_from_path(self, file_path):
        image = io.imread(file_path)

        self.image_min = np.percentile(image, 1)
        self.image_max = np.percentile(image, 99)
        self.resolution = (self.image_max - self.image_min) / 100 if self.image_max != self.image_min else 1

        self.intensity_threshold_slider.config(
            from_=self.image_min,
            to=self.image_max,
            resolution=self.resolution,
        )

        self.file_path = file_path
        self.show_image(file_path)

        self.checkbox.config(state="normal")
        self.save_button.config(state="normal")

        if self.show_size_sliders:
            self.min_size_threshold_slider.config(state="normal")
            self.max_size_threshold_slider.config(state="normal")
            self.min_size_threshold_slider.config(troughcolor=self.troughcolor)
            self.max_size_threshold_slider.config(troughcolor=self.troughcolor)
            self.min_size_threshold_slider.config(fg=self.fg)
            self.max_size_threshold_slider.config(fg=self.fg)

        self.checkbox.config(fg=self.fg)

        try:
            self.intensity_threshold = int(self.intensity_threshold)
        except ValueError:
            pass

        if isinstance(self.intensity_threshold, (int, float)):
            self.intensity_threshold_slider.config(state="normal")
            self.intensity_threshold_slider.set(self.intensity_threshold)
            self.intensity_threshold_slider.config(state="disabled")

        self.apply_checkbox()

    def show_image(self, file_path):
        if self.show_size_sliders:
            min_size_pixels = self.min_size_threshold * (self.scale ** 2)
            max_size_pixels = self.max_size_threshold * (self.scale ** 2)
            image, _ = intensity_threshold(
                file_path,
                min_size_pixels,
                max_size_pixels,
                self.intensity_threshold,
            )
        else:
            # Cellbody preview: intensity threshold only, no watershed and no size filtering.
            image, _, _ = threshold_cellbody_without_size_filter(
                file_path,
                self.intensity_threshold,
            )

        self.image = image
        self.set_preview_image(image)


class CellSizePreviewWindow(BasePreviewWindow):
    def __init__(
        self,
        master,
        scale,
        intensity_threshold_cell,
        default_min,
        default_max,
        watershed_mode="No Watershed",
        watershed_sigma=6.0,
        watershed_min_distance=30,
        watershed_h_value=3.0,
        initial_file_path=None,
        cached_raw_stage=None,
    ):
        self.master = master
        self.setup_window("Cell Size Preview", height=720)

        self.scale = scale
        self.intensity_threshold_cell = intensity_threshold_cell
        self.min_size_threshold = default_min
        self.max_size_threshold = default_max
        self.watershed_mode = watershed_mode
        self.watershed_sigma = watershed_sigma
        self.watershed_min_distance = watershed_min_distance
        self.watershed_h_value = watershed_h_value
        self.file_path = initial_file_path
        self.cached_raw_stage = cached_raw_stage
        self._raw_stage_cache_key = None
        self._raw_labels_cache = None
        self._raw_info_cache = None

        self.image_label = self.create_image_label()

        info = "Preview shows the cellbody objects after cell size and border filtering."
        self.info_label = ttk.Label(self.master, text=info, wraplength=self.WINDOW_WIDTH - 25, justify="left")
        self.info_label.pack(pady=(0, 4), padx=10)

        self.load_button = ttk.Button(self.master, text="Load Image", command=self.load_image)
        self.load_button.pack(pady=5)

        self.min_size_threshold_slider = self.create_slider(
            from_=0,
            to=54,
            resolution=1,
            command=self.update_parameter,
        )
        self.min_size_threshold_slider.set(slider_position_value(self.min_size_threshold))

        self.max_size_threshold_slider = self.create_slider(
            from_=0,
            to=54,
            resolution=1,
            command=self.update_parameter,
        )
        self.max_size_threshold_slider.set(slider_position_value(self.max_size_threshold))

        self.min_size_threshold_slider.config(showvalue=0)
        self.max_size_threshold_slider.config(showvalue=0)

        self.save_button = ttk.Button(self.master, text="Save Size Thresholds", command=self.save_threshold)
        self.save_button.pack(pady=5)

        self.create_blank_preview()
        self.update_parameter(0)

        if self.file_path:
            self.show_image(self.file_path)

    def save_threshold(self):
        self.master.destroy()

    def load_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.file_path = file_path
            self.show_image(file_path)

    def update_parameter(self, _value):
        self.min_size_threshold = stepped_slider_sequence(self.min_size_threshold_slider.get())
        self.max_size_threshold = stepped_slider_sequence(self.max_size_threshold_slider.get())

        self.min_size_threshold_slider.config(label=f"Min Size Threshold: {self.min_size_threshold:.0f} µm²")
        self.max_size_threshold_slider.config(label=f"Max Size Threshold: {self.max_size_threshold:.0f} µm²")

        if self.file_path:
            self.show_image(self.file_path)

    def _make_raw_stage_key(self, file_path):
        return (
            os.path.abspath(file_path) if file_path else None,
            str(self.intensity_threshold_cell),
            self.watershed_mode,
            round(float(self.watershed_sigma), 4),
            int(self.watershed_min_distance),
            round(float(self.watershed_h_value), 4),
        )

    def _copy_info(self, info):
        return {int(k): dict(v) for k, v in (info or {}).items()}

    def _prepare_raw_stage(self, file_path):
        key = self._make_raw_stage_key(file_path)

        if self._raw_stage_cache_key == key and self._raw_labels_cache is not None:
            return

        if (
            isinstance(self.cached_raw_stage, dict)
            and self.cached_raw_stage.get("key") == key
            and self.cached_raw_stage.get("raw_labels") is not None
        ):
            self._raw_labels_cache = np.array(self.cached_raw_stage["raw_labels"], copy=True)
            self._raw_info_cache = self._copy_info(self.cached_raw_stage.get("raw_info", {}))
            self._raw_stage_cache_key = key
            return

        _, _, binary_mask_filled = threshold_cellbody_without_size_filter(
            file_path,
            self.intensity_threshold_cell,
        )

        if self.watershed_mode == "No Watershed":
            raw_labels = measure.label(binary_mask_filled)
            raw_info = {}
        else:
            raw_labels, raw_info = build_watershed_labels_without_size_filter(
                binary_mask_filled,
                sigma=self.watershed_sigma,
                min_distance=self.watershed_min_distance,
                h_value=self.watershed_h_value,
            )

        self._raw_labels_cache = raw_labels
        self._raw_info_cache = self._copy_info(raw_info)
        self._raw_stage_cache_key = key

    def show_image(self, file_path):
        min_size_pixels = self.min_size_threshold * (self.scale ** 2)
        max_size_pixels = self.max_size_threshold * (self.scale ** 2)

        self._prepare_raw_stage(file_path)
        raw_labels = self._raw_labels_cache
        raw_info = self._raw_info_cache

        if self.watershed_mode == "Exclude":
            not_split_labels, not_split_info = subset_label_image_by_predicate(
                raw_labels,
                raw_info,
                lambda info: int(info.get("raw_number_of_segments_in_original_object", 1)) <= 1,
            )
            analysis_labels, analysis_info, _removed = size_filter_label_image(
                not_split_labels,
                not_split_info,
                min_size_pixels,
                max_size_pixels,
            )
        else:
            analysis_labels, analysis_info, _removed = size_filter_label_image(
                raw_labels,
                raw_info,
                min_size_pixels,
                max_size_pixels,
            )

        self.image = labels_to_gray_display(analysis_labels, analysis_info)
        self.set_preview_image(self.image)

class WatershedPreviewWindow(BasePreviewWindow):
    def __init__(
        self,
        master,
        scale,
        intensity_threshold_cell,
        min_size_cell,
        max_size_cell,
        default_sigma,
        default_min_distance,
        default_h_value,
        watershed_mode="Split",
        initial_file_path=None,
    ):
        self.master = master
        self.setup_window("Watershed Preview", height=730)

        self.scale = scale
        self.intensity_threshold_cell = intensity_threshold_cell
        self.min_size_cell = min_size_cell
        self.max_size_cell = max_size_cell
        self.watershed_mode = watershed_mode

        self.sigma = float(default_sigma)
        self.min_distance = int(default_min_distance)
        self.h_value = float(default_h_value)

        self.file_path = initial_file_path
        self.last_watershed_cache = None

        self._preview_after_id = None
        self._is_rendering = False

        self.image_label = self.create_image_label()

        info_text = (
            "Split: preview shows watershed-separated cellbody objects before size filtering."
            if watershed_mode == "Split"
            else "Exclude: preview shows watershed-affected cellbody objects before exclusion and size filtering."
        )
        self.info_label = ttk.Label(self.master, text=info_text, wraplength=self.WINDOW_WIDTH - 25, justify="left")
        self.info_label.pack(pady=(0, 4), padx=10)

        self.load_button = ttk.Button(self.master, text="Load Image", command=self.load_image)
        self.load_button.pack(pady=5)

        self.sigma_slider = self.create_slider(
            from_=0,
            to=20,
            resolution=0.1,
            command=self.update_parameter,
            label=f"Sigma: {self.sigma:.1f}",
        )
        self.sigma_slider.set(self.sigma)

        self.min_distance_slider = self.create_slider(
            from_=1,
            to=100,
            resolution=1,
            command=self.update_parameter,
            label=f"Min Distance: {self.min_distance}",
        )
        self.min_distance_slider.set(self.min_distance)

        self.h_value_slider = self.create_slider(
            from_=0,
            to=20,
            resolution=0.1,
            command=self.update_parameter,
            label=f"h-Value: {self.h_value:.1f}",
        )
        self.h_value_slider.set(self.h_value)

        self.save_button = ttk.Button(
            self.master,
            text="Save Watershed Parameters",
            command=self.save_parameters,
        )
        self.save_button.pack(pady=4)

        self.create_blank_preview()

        if self.file_path:
            self.schedule_preview_update()

    def save_parameters(self):
        self.master.destroy()

    def load_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.file_path = file_path
            self.schedule_preview_update()

    def update_parameter(self, _value=None):
        self.sigma = float(self.sigma_slider.get())
        self.min_distance = int(self.min_distance_slider.get())
        self.h_value = float(self.h_value_slider.get())

        self.sigma_slider.config(label=f"Sigma: {self.sigma:.1f}")
        self.min_distance_slider.config(label=f"Min Distance: {self.min_distance}")
        self.h_value_slider.config(label=f"h-Value: {self.h_value:.1f}")

        self.schedule_preview_update()

    def schedule_preview_update(self):
        if not self.file_path:
            return

        if self._preview_after_id is not None:
            self.master.after_cancel(self._preview_after_id)

        self._preview_after_id = self.master.after(90, self.show_image)

    def show_image(self):
        if not self.file_path or self._is_rendering:
            return

        self._preview_after_id = None
        self._is_rendering = True

        try:
            _, _, binary_mask_filled = threshold_cellbody_without_size_filter(
                self.file_path,
                self.intensity_threshold_cell,
            )
            raw_labels, raw_info = build_watershed_labels_without_size_filter(
                binary_mask_filled,
                sigma=self.sigma,
                min_distance=self.min_distance,
                h_value=self.h_value,
            )
            cache_key = (
                os.path.abspath(self.file_path) if self.file_path else None,
                str(self.intensity_threshold_cell),
                self.watershed_mode,
                round(float(self.sigma), 4),
                int(self.min_distance),
                round(float(self.h_value), 4),
            )
            self.last_watershed_cache = {
                "key": cache_key,
                "raw_labels": np.array(raw_labels, copy=True),
                "raw_info": {int(k): dict(v) for k, v in raw_info.items()},
            }

            preview = labels_to_gray_display(raw_labels, raw_info)
            self.set_preview_image(preview)

        except Exception as exc:
            messagebox.showerror("Watershed Preview Error", str(exc))
        finally:
            self._is_rendering = False




class RatioPreviewWindow(BasePreviewWindow):
    def __init__(
        self,
        master,
        scale,
        intensity_threshold_cell,
        min_size_cell,
        max_size_cell,
        intensity_threshold_nucleus,
        min_size_nucleus,
        max_size_nucleus,
        circularity_threshold,
        watershed_mode,
        watershed_sigma,
        watershed_min_distance,
        watershed_h_value,
        ratio_threshold,
        initial_cell_file_path=None,
        initial_nucleus_file_path=None,
    ):
        self.master = master
        self.setup_window("Ratio Preview", height=580)

        self.scale = scale
        self.intensity_threshold_cell = intensity_threshold_cell
        self.min_size_cell = min_size_cell
        self.max_size_cell = max_size_cell
        self.intensity_threshold_nucleus = intensity_threshold_nucleus
        self.min_size_nucleus = min_size_nucleus
        self.max_size_nucleus = max_size_nucleus
        self.circularity_threshold = circularity_threshold

        self.watershed_mode = watershed_mode
        self.watershed_sigma = watershed_sigma
        self.watershed_min_distance = watershed_min_distance
        self.watershed_h_value = watershed_h_value

        self.ratio_threshold = float(ratio_threshold)

        self.cell_file_path = initial_cell_file_path
        self.nucleus_file_path = initial_nucleus_file_path

        self._cached_valid_objects = None
        self._cached_preview_shape = None
        self._cached_label_image = None
        self._cached_watershed_segment_info = None
        self._cache_key = None
        self._preview_after_id = None
        self._is_rendering = False

        self.image_label = self.create_image_label()

        button_frame = ttk.Frame(self.master)
        button_frame.pack(pady=5)

        self.load_cell_button = ttk.Button(
            button_frame,
            text="Load Cellbody Image",
            command=self.load_cell_image,
        )
        self.load_cell_button.pack(side="left", padx=5)

        self.load_nucleus_button = ttk.Button(
            button_frame,
            text="Load Nucleus Image",
            command=self.load_nucleus_image,
        )
        self.load_nucleus_button.pack(side="left", padx=5)

        self.ratio_slider = self.create_slider(
            from_=0,
            to=1,
            resolution=0.01,
            command=self.update_parameter,
            label=f"Nucleus / Cellbody Area Ratio < {self.ratio_threshold:.2f}",
        )
        self.ratio_slider.set(self.ratio_threshold)

        self.save_button = ttk.Button(
            self.master,
            text="Save Ratio Threshold",
            command=self.save_parameters,
        )
        self.save_button.pack(pady=8)

        self.create_blank_preview()

        if self.cell_file_path and self.nucleus_file_path:
            self.schedule_preview_update(rebuild_cache=True)

    def save_parameters(self):
        self.master.destroy()

    def load_cell_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.cell_file_path = file_path
            self.schedule_preview_update(rebuild_cache=True)

    def load_nucleus_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.nucleus_file_path = file_path
            self.schedule_preview_update(rebuild_cache=True)

    def update_parameter(self, _value=None):
        self.ratio_threshold = float(self.ratio_slider.get())
        self.ratio_slider.config(
            label=f"Nucleus / Cellbody Area Ratio < {self.ratio_threshold:.2f}"
        )
        self.schedule_preview_update(rebuild_cache=False)

    def build_cache_key(self):
        return (
            self.cell_file_path,
            self.nucleus_file_path,
            self.scale,
            self.intensity_threshold_cell,
            self.min_size_cell,
            self.max_size_cell,
            self.intensity_threshold_nucleus,
            self.min_size_nucleus,
            self.max_size_nucleus,
            self.circularity_threshold,
            self.watershed_mode,
            self.watershed_sigma,
            self.watershed_min_distance,
            self.watershed_h_value,
        )

    def ensure_cached_objects(self):
        cache_key = self.build_cache_key()

        if self._cache_key == cache_key and self._cached_valid_objects is not None:
            return

        valid_objects, preview_shape, label_image, watershed_segment_info = ratio_filter_preview(
            self.cell_file_path,
            self.nucleus_file_path,
            self.scale,
            self.intensity_threshold_cell,
            self.min_size_cell,
            self.max_size_cell,
            self.intensity_threshold_nucleus,
            self.min_size_nucleus,
            self.max_size_nucleus,
            self.circularity_threshold,
            self.watershed_mode,
            self.watershed_sigma,
            self.watershed_min_distance,
            self.watershed_h_value,
        )

        self._cached_valid_objects = valid_objects
        self._cached_preview_shape = preview_shape
        self._cached_label_image = label_image
        self._cached_watershed_segment_info = watershed_segment_info
        self._cache_key = cache_key

    def schedule_preview_update(self, rebuild_cache=False):
        if not self.cell_file_path or not self.nucleus_file_path:
            return

        if rebuild_cache:
            self._cached_valid_objects = None
            self._cached_preview_shape = None
            self._cached_label_image = None
            self._cached_watershed_segment_info = None
            self._cache_key = None

        if self._preview_after_id is not None:
            self.master.after_cancel(self._preview_after_id)

        self._preview_after_id = self.master.after(90, self.show_image)

    def show_image(self):
        if not self.cell_file_path or not self.nucleus_file_path or self._is_rendering:
            return

        self._preview_after_id = None
        self._is_rendering = True

        try:
            self.ensure_cached_objects()

            keep_labels = set()
            nucleus_holes = np.zeros(self._cached_preview_shape, dtype=bool)

            for obj in self._cached_valid_objects:
                cell_area = obj["cell_area"]
                nucleus_area = obj["nucleus_area"]

                if cell_area <= 0:
                    continue

                if nucleus_area / cell_area < self.ratio_threshold:
                    keep_labels.add(int(obj["label"]))
                    nucleus_coords = obj["nucleus_coords"]
                    nucleus_holes[nucleus_coords[:, 0], nucleus_coords[:, 1]] = True

            ratio_label_image = np.zeros_like(self._cached_label_image, dtype=int)
            for label in keep_labels:
                ratio_label_image[self._cached_label_image == label] = label

            preview = labels_to_gray_display(
                ratio_label_image,
                self._cached_watershed_segment_info,
                nucleus_holes_mask=nucleus_holes,
            )

            self.set_preview_image(preview)

        except Exception as exc:
            messagebox.showerror("Ratio Preview Error", str(exc))
        finally:
            self._is_rendering = False



class CircularityPreviewWindow(BasePreviewWindow):
    def __init__(self, master, nucleus_image):
        self.nucleus_image = nucleus_image
        self.circularity_threshold = 0.73
        self.master = master
        self.setup_window("Image Viewer", height=540)

        self.image_label = self.create_image_label()

        self.circularity_threshold_slider = self.create_slider(
            from_=0,
            to=1,
            resolution=0.01,
            command=self.update_parameter,
        )
        self.circularity_threshold_slider.set(0.73)

        self.save_button = ttk.Button(self.master, text="Save Threshold", command=self.save_threshold)
        self.save_button.pack(pady=5)

        self.set_preview_image(self.nucleus_image)

    def save_threshold(self):
        self.master.destroy()

    def update_parameter(self, _value):
        self.circularity_threshold = self.circularity_threshold_slider.get()
        self.circularity_threshold_slider.config(
            label=f"Circularity Threshold: {self.circularity_threshold:.2f}"
        )
        self.show_image()

    def show_image(self):
        circularity_mask = circularity_filter_preview(self.nucleus_image, self.circularity_threshold)
        self.set_preview_image(circularity_mask)


# Parameter GUI

class ParameterGUI:
    def __init__(self, defaults):
        self.defaults = defaults.copy()
        self.parameters = None
        self.nucleus_image = None
        self.cell_preview_image_path = None
        self.nucleus_preview_image_path = None
        self.cellbody_raw_stage_cache = None

        self.watershed_sigma = self.defaults["watershed_sigma"]
        self.watershed_min_distance = self.defaults["watershed_min_distance"]
        self.watershed_h_value = self.defaults["watershed_h_value"]

        self.root = tk.Tk()
        self.root.title("ViSCA - Image Setting Parameters")
        self.root.geometry("500x760")
        self.root.minsize(465, 650)

        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True)

        self.canvas = tk.Canvas(outer, borderwidth=0, highlightthickness=0)
        scrollbar = ttk.Scrollbar(outer, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.parameters_frame = ttk.Frame(self.canvas, padding="10")
        self.canvas_window = self.canvas.create_window((0, 0), window=self.parameters_frame, anchor="nw")

        def _on_frame_configure(_event):
            self.canvas.configure(scrollregion=self.canvas.bbox("all"))

        def _on_canvas_configure(event):
            self.canvas.itemconfigure(self.canvas_window, width=event.width)

        def _on_mousewheel(event):
            if getattr(event, "num", None) == 4:
                self.canvas.yview_scroll(-4, "units")
            elif getattr(event, "num", None) == 5:
                self.canvas.yview_scroll(4, "units")
            else:
                delta = getattr(event, "delta", 0)
                if delta == 0:
                    return "break"
                if abs(delta) < 120:
                    step = -4 if delta > 0 else 4
                else:
                    step = int(-delta / 120) * 4
                self.canvas.yview_scroll(step, "units")
            return "break"

        def _start_drag_scroll(event):
            self.canvas.scan_mark(event.x_root, event.y_root)

        def _drag_scroll(event):
            self.canvas.scan_dragto(event.x_root, event.y_root, gain=1)
            return "break"

        def _bind_scroll_recursive(widget):
            widget.bind("<MouseWheel>", _on_mousewheel, add="+")
            widget.bind("<Button-4>", _on_mousewheel, add="+")
            widget.bind("<Button-5>", _on_mousewheel, add="+")
            for child in widget.winfo_children():
                _bind_scroll_recursive(child)

        self.parameters_frame.bind("<Configure>", _on_frame_configure)
        self.canvas.bind("<Configure>", _on_canvas_configure)
        self.canvas.bind("<MouseWheel>", _on_mousewheel, add="+")
        self.canvas.bind("<Button-4>", _on_mousewheel, add="+")
        self.canvas.bind("<Button-5>", _on_mousewheel, add="+")
        self.parameters_frame.bind("<MouseWheel>", _on_mousewheel, add="+")
        self.parameters_frame.bind("<Button-4>", _on_mousewheel, add="+")
        self.parameters_frame.bind("<Button-5>", _on_mousewheel, add="+")
        self.canvas.bind("<ButtonPress-1>", _start_drag_scroll, add="+")
        self.canvas.bind("<B1-Motion>", _drag_scroll, add="+")
        self.root.bind_all("<MouseWheel>", _on_mousewheel, add="+")
        self.root.bind_all("<Button-4>", _on_mousewheel, add="+")
        self.root.bind_all("<Button-5>", _on_mousewheel, add="+")
        self.root.after(200, lambda: _bind_scroll_recursive(self.parameters_frame))

        bottom_frame = ttk.Frame(self.root, padding=(10, 6))
        bottom_frame.pack(fill="x")
        ttk.Label(
            bottom_frame,
            text="Scroll down for all parameters ↓",
            foreground="dark grey",
        ).pack(anchor="w", pady=(0, 4))
        self.save_button = ttk.Button(
            bottom_frame,
            text="Save Parameters & Start Run",
            command=self.save_parameters,
        )
        self.save_button.pack(fill="x")

        row = 0
        title_label = ttk.Label(
            self.parameters_frame,
            text="ViSCA",
            font=("Arial", 20, "bold"),
            foreground=required_color,
            anchor="center",
        )
        title_label.grid(row=row, column=0, columnspan=2, sticky="ew")
        row += 1

        main_heading = ttk.Label(
            self.parameters_frame,
            text="Set Parameters or use default settings",
            font=("Arial", 16, "bold"),
        )
        main_heading.grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(8, 4))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="Scroll down for all parameters ↓",
            foreground="dark grey",
            wraplength=590,
            justify="left",
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 8))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="*Scale (pixels per micron):",
            foreground=required_color,
        ).grid(row=row, column=0, sticky=tk.W)
        self.scale_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.scale_entry, str(self.defaults["scale"]))
        self.scale_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="Cellbody Parameters",
            font=("Arial", 13, "bold"),
        ).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10, 3))
        row += 1

        ttk.Button(
            self.parameters_frame,
            text="Adjust intensity threshold using one of your cellbody images",
            command=self.open_parameter_preview_cell,
        ).grid(row=row, column=0, columnspan=2, sticky="", pady=(2, 2))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="*Intensity Threshold cellbody:",
            foreground=required_color,
        ).grid(row=row, column=0, sticky=tk.W)
        self.intensity_threshold_cell_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.intensity_threshold_cell_entry, str(self.defaults["intensity_threshold_cell"]))
        self.intensity_threshold_cell_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1


        ttk.Label(
            self.parameters_frame,
            text="*Watershed mode:",
            foreground=required_color,
        ).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10, 3))
        row += 1

        self.watershed_mode_var = tk.StringVar(value=self.defaults.get("watershed_mode", "Split"))
        mode_frame = ttk.Frame(self.parameters_frame)
        mode_frame.grid(row=row, column=0, columnspan=2, sticky="w")
        for mode in ["No Watershed", "Split", "Exclude"]:
            ttk.Radiobutton(
                mode_frame,
                text=mode,
                value=mode,
                variable=self.watershed_mode_var,
                command=self.toggle_watershed_widgets,
            ).pack(side="left", padx=(0, 14))
        row += 1

        self.watershed_static_info_label = ttk.Label(
            self.parameters_frame,
            text=(
                "Split: split connected cell-body objects and keep valid segments\n"
                "Exclude: detect watershed-affected objects and exclude them"
            ),
            foreground="dark grey",
            wraplength=590,
            justify="left",
        )
        self.watershed_static_info_label.grid(row=row, column=0, columnspan=2, sticky="w", pady=(2, 4))
        row += 1

        self.watershed_button = ttk.Button(
            self.parameters_frame,
            text="Adjust watershed using one of your cellbody images",
            command=self.open_watershed_preview,
        )
        self.watershed_button_row = row
        row += 1

        self.watershed_sigma_label = ttk.Label(self.parameters_frame, text="Watershed sigma:")
        self.watershed_sigma_entry = ttk.Entry(self.parameters_frame)
        self.watershed_sigma_entry.insert(0, str(self.watershed_sigma))
        self.watershed_sigma_row = row
        row += 1

        self.watershed_min_distance_label = ttk.Label(self.parameters_frame, text="Watershed min distance:")
        self.watershed_min_distance_entry = ttk.Entry(self.parameters_frame)
        self.watershed_min_distance_entry.insert(0, str(self.watershed_min_distance))
        self.watershed_min_distance_row = row
        row += 1

        self.watershed_h_value_label = ttk.Label(self.parameters_frame, text="Watershed h-value:")
        self.watershed_h_value_entry = ttk.Entry(self.parameters_frame)
        self.watershed_h_value_entry.insert(0, str(self.watershed_h_value))
        self.watershed_h_value_row = row
        row += 1

        ttk.Button(
            self.parameters_frame,
            text="Adjust cell size filter using one of your cellbody images",
            command=self.open_cell_size_preview,
        ).grid(row=row, column=0, columnspan=2, sticky="", pady=(4, 0))
        row += 1

        ttk.Label(self.parameters_frame, text="Min size cell [µm²]:").grid(row=row, column=0, sticky=tk.W)
        self.min_size_cell_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.min_size_cell_entry, str(self.defaults["min_size_cell"]))
        self.min_size_cell_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Label(self.parameters_frame, text="Max size cell [µm²]:").grid(row=row, column=0, sticky=tk.W)
        self.max_size_cell_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.max_size_cell_entry, str(self.defaults["max_size_cell"]))
        self.max_size_cell_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="Nucleus Parameters",
            font=("Arial", 13, "bold"),
        ).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10, 3))
        row += 1

        ttk.Button(
            self.parameters_frame,
            text="Adjust intensity threshold using one of your nucleus images",
            command=self.open_parameter_preview_nucleus,
        ).grid(row=row, column=0, columnspan=2, sticky="", pady=(2, 2))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="*Intensity Threshold Nucleus:",
            foreground=required_color,
        ).grid(row=row, column=0, sticky=tk.W)
        self.intensity_threshold_nucleus_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.intensity_threshold_nucleus_entry, str(self.defaults["intensity_threshold_nucleus"]))
        self.intensity_threshold_nucleus_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Label(self.parameters_frame, text="Min size nucleus [µm²]:").grid(row=row, column=0, sticky=tk.W)
        self.min_size_nucleus_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.min_size_nucleus_entry, str(self.defaults["min_size_nucleus"]))
        self.min_size_nucleus_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Label(self.parameters_frame, text="Max size nucleus [µm²]:").grid(row=row, column=0, sticky=tk.W)
        self.max_size_nucleus_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.max_size_nucleus_entry, str(self.defaults["max_size_nucleus"]))
        self.max_size_nucleus_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Button(
            self.parameters_frame,
            text="Adjust circularity threshold using one of your nucleus images",
            command=self.open_parameter_preview_circularity,
        ).grid(row=row, column=0, columnspan=2, sticky="", pady=(10, 0))
        row += 1

        ttk.Label(self.parameters_frame, text="Exclude cells with nucleus circularity ≤:").grid(row=row, column=0, sticky=tk.W)
        self.circularity_threshold_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.circularity_threshold_entry, str(self.defaults["circularity_threshold"]))
        self.circularity_threshold_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Button(
            self.parameters_frame,
            text="Adjust nucleus/cellbody area ratio using your example images",
            command=self.open_ratio_preview,
        ).grid(row=row, column=0, columnspan=2, sticky="", pady=(10, 0))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="Exclude cells with nucleus/cellbody area ratio ≥:",
        ).grid(row=row, column=0, sticky=tk.W)
        self.ratio_entry = ttk.Entry(self.parameters_frame)
        add_placeholder(self.ratio_entry, str(self.defaults["ratio_threshold"]))
        self.ratio_entry.grid(row=row, column=1, sticky=(tk.W, tk.E))
        row += 1

        ttk.Label(
            self.parameters_frame,
            text="* = This parameter has to be adjusted. Don't rely on the default value!",
            foreground=required_color,
        ).grid(row=row, column=0, columnspan=2, sticky="ew", pady=(10, 3))

        self.toggle_watershed_widgets()

        self.parameters_frame.columnconfigure(1, weight=1)

    def toggle_watershed_widgets(self):
        if self.watershed_mode_var.get() != "No Watershed":
            self.watershed_button.grid(row=self.watershed_button_row, column=0, columnspan=2, sticky="", pady=(2, 0))
            self.watershed_sigma_label.grid(row=self.watershed_sigma_row, column=0, sticky=tk.W)
            self.watershed_sigma_entry.grid(row=self.watershed_sigma_row, column=1, sticky=(tk.W, tk.E))
            self.watershed_min_distance_label.grid(row=self.watershed_min_distance_row, column=0, sticky=tk.W)
            self.watershed_min_distance_entry.grid(row=self.watershed_min_distance_row, column=1, sticky=(tk.W, tk.E))
            self.watershed_h_value_label.grid(row=self.watershed_h_value_row, column=0, sticky=tk.W)
            self.watershed_h_value_entry.grid(row=self.watershed_h_value_row, column=1, sticky=(tk.W, tk.E))
        else:
            self.watershed_button.grid_remove()
            self.watershed_sigma_label.grid_remove()
            self.watershed_sigma_entry.grid_remove()
            self.watershed_min_distance_label.grid_remove()
            self.watershed_min_distance_entry.grid_remove()
            self.watershed_h_value_label.grid_remove()
            self.watershed_h_value_entry.grid_remove()

    def _get_float(self, entry):
        return float(entry.get())

    def _get_int(self, entry):
        return int(entry.get())

    def read_parameters_from_entries(self):
        watershed_mode = self.watershed_mode_var.get()
        watershed_sigma = self.watershed_sigma
        watershed_min_distance = self.watershed_min_distance
        watershed_h_value = self.watershed_h_value

        if watershed_mode != "No Watershed":
            watershed_sigma = float(self.watershed_sigma_entry.get())
            watershed_min_distance = int(self.watershed_min_distance_entry.get())
            watershed_h_value = float(self.watershed_h_value_entry.get())

        return {
            "scale": self._get_float(self.scale_entry),
            "intensity_threshold_cell": self.intensity_threshold_cell_entry.get(),
            "min_size_cell": self._get_int(self.min_size_cell_entry),
            "max_size_cell": self._get_int(self.max_size_cell_entry),
            "intensity_threshold_nucleus": self.intensity_threshold_nucleus_entry.get(),
            "min_size_nucleus": self._get_int(self.min_size_nucleus_entry),
            "max_size_nucleus": self._get_int(self.max_size_nucleus_entry),
            "circularity_threshold": self._get_float(self.circularity_threshold_entry),
            "ratio_threshold": self._get_float(self.ratio_entry),
            "watershed_mode": watershed_mode,
            "watershed_sigma": watershed_sigma,
            "watershed_min_distance": watershed_min_distance,
            "watershed_h_value": watershed_h_value,
        }

    def save_parameters(self):
        try:
            parameters = self.read_parameters_from_entries()
            validate_parameters(parameters)
            self.parameters = parameters

            messagebox.showinfo(
                "Parameters Saved",
                "Parameters have been saved. Now select images of cell bodies.",
            )
            self.root.destroy()

        except ValueError:
            messagebox.showerror(
                "Input Error",
                "Please ensure all fields are filled correctly with numerical values.",
            )

    def open_parameter_preview_cell(self):
        try:
            scale = float(self.scale_entry.get())
            intensity_threshold_cell = self.intensity_threshold_cell_entry.get()
            min_size_cell = int(self.min_size_cell_entry.get())
            max_size_cell = int(self.max_size_cell_entry.get())
        except ValueError:
            messagebox.showerror("Input Error", "Please fill scale and cellbody threshold/size parameters first.")
            return

        parameter_preview_window = tk.Toplevel(self.root)
        parameter_preview = ParameterPreviewWindow(
            parameter_preview_window,
            scale,
            intensity_threshold_cell,
            min_size_cell,
            max_size_cell,
            "mean",
            initial_file_path=self.cell_preview_image_path,
            show_size_sliders=False,
        )
        self.root.wait_window(parameter_preview_window)

        set_entry_value(self.intensity_threshold_cell_entry, parameter_preview.intensity_threshold)
        self.cell_preview_image_path = parameter_preview.file_path
        self.cellbody_raw_stage_cache = None

    def open_cell_size_preview(self):
        try:
            scale = float(self.scale_entry.get())
            intensity_threshold_cell = self.intensity_threshold_cell_entry.get()
            min_size_cell = int(self.min_size_cell_entry.get())
            max_size_cell = int(self.max_size_cell_entry.get())
            watershed_mode = self.watershed_mode_var.get()

            if watershed_mode != "No Watershed":
                watershed_sigma = float(self.watershed_sigma_entry.get())
                watershed_min_distance = int(self.watershed_min_distance_entry.get())
                watershed_h_value = float(self.watershed_h_value_entry.get())
            else:
                watershed_sigma = self.watershed_sigma
                watershed_min_distance = self.watershed_min_distance
                watershed_h_value = self.watershed_h_value
        except ValueError:
            messagebox.showerror("Input Error", "Please fill scale, cellbody, and watershed parameters first.")
            return

        size_preview_window = tk.Toplevel(self.root)
        size_preview = CellSizePreviewWindow(
            size_preview_window,
            scale,
            intensity_threshold_cell,
            min_size_cell,
            max_size_cell,
            watershed_mode=watershed_mode,
            watershed_sigma=watershed_sigma,
            watershed_min_distance=watershed_min_distance,
            watershed_h_value=watershed_h_value,
            initial_file_path=self.cell_preview_image_path,
            cached_raw_stage=self.cellbody_raw_stage_cache,
        )
        self.root.wait_window(size_preview_window)

        set_entry_value(self.min_size_cell_entry, int(size_preview.min_size_threshold))
        set_entry_value(self.max_size_cell_entry, int(size_preview.max_size_threshold))
        self.cell_preview_image_path = size_preview.file_path
        self.cellbody_raw_stage_cache = {
            "key": size_preview._raw_stage_cache_key,
            "raw_labels": np.array(size_preview._raw_labels_cache, copy=True) if size_preview._raw_labels_cache is not None else None,
            "raw_info": {int(k): dict(v) for k, v in (size_preview._raw_info_cache or {}).items()},
        } if size_preview._raw_stage_cache_key is not None else self.cellbody_raw_stage_cache

    def open_watershed_preview(self):
        if self.watershed_mode_var.get() == "No Watershed":
            return

        try:
            scale = float(self.scale_entry.get())
            intensity_threshold_cell = self.intensity_threshold_cell_entry.get()
            sigma = float(self.watershed_sigma_entry.get())
            min_distance = int(self.watershed_min_distance_entry.get())
            h_value = float(self.watershed_h_value_entry.get())
        except ValueError:
            messagebox.showerror(
                "Input Error",
                "Please fill scale, cellbody threshold, and watershed parameters first."
            )
            return

        watershed_preview_window = tk.Toplevel(self.root)
        watershed_preview = WatershedPreviewWindow(
            watershed_preview_window,
            scale,
            intensity_threshold_cell,
            0,
            0,
            sigma,
            min_distance,
            h_value,
            watershed_mode=self.watershed_mode_var.get(),
            initial_file_path=self.cell_preview_image_path,
        )
        self.root.wait_window(watershed_preview_window)

        self.watershed_sigma = watershed_preview.sigma
        self.watershed_min_distance = watershed_preview.min_distance
        self.watershed_h_value = watershed_preview.h_value

        set_entry_value(self.watershed_sigma_entry, self.watershed_sigma)
        set_entry_value(self.watershed_min_distance_entry, self.watershed_min_distance)
        set_entry_value(self.watershed_h_value_entry, self.watershed_h_value)

        if watershed_preview.file_path:
            self.cell_preview_image_path = watershed_preview.file_path

        if getattr(watershed_preview, "last_watershed_cache", None) is not None:
            self.cellbody_raw_stage_cache = watershed_preview.last_watershed_cache

    def open_parameter_preview_nucleus(self):
        try:
            scale = float(self.scale_entry.get())
            intensity_threshold_nucleus = self.intensity_threshold_nucleus_entry.get()
            min_size_nucleus = int(self.min_size_nucleus_entry.get())
            max_size_nucleus = int(self.max_size_nucleus_entry.get())
        except ValueError:
            messagebox.showerror("Input Error", "Please fill scale and nucleus parameters first.")
            return

        parameter_preview_window = tk.Toplevel(self.root)
        parameter_preview = ParameterPreviewWindow(
            parameter_preview_window,
            scale,
            intensity_threshold_nucleus,
            min_size_nucleus,
            max_size_nucleus,
            "otsu",
            initial_file_path=self.nucleus_preview_image_path,
        )
        self.root.wait_window(parameter_preview_window)

        self.nucleus_image = parameter_preview.image

        set_entry_value(self.intensity_threshold_nucleus_entry, parameter_preview.intensity_threshold)
        set_entry_value(self.min_size_nucleus_entry, int(parameter_preview.min_size_threshold))
        set_entry_value(self.max_size_nucleus_entry, int(parameter_preview.max_size_threshold))

        self.nucleus_preview_image_path = parameter_preview.file_path

    def open_parameter_preview_circularity(self):
        if self.nucleus_image is None:
            messagebox.showerror(
                "Input Error",
                "Please adjust nucleus parameters first so a nucleus preview image is available.",
            )
            return

        circularity_preview_window = tk.Toplevel(self.root)
        circularity_preview = CircularityPreviewWindow(circularity_preview_window, self.nucleus_image)
        self.root.wait_window(circularity_preview_window)

        set_entry_value(self.circularity_threshold_entry, circularity_preview.circularity_threshold)

    def open_ratio_preview(self):
        try:
            scale = float(self.scale_entry.get())
            intensity_threshold_cell = self.intensity_threshold_cell_entry.get()
            min_size_cell = int(self.min_size_cell_entry.get())
            max_size_cell = int(self.max_size_cell_entry.get())

            intensity_threshold_nucleus = self.intensity_threshold_nucleus_entry.get()
            min_size_nucleus = int(self.min_size_nucleus_entry.get())
            max_size_nucleus = int(self.max_size_nucleus_entry.get())

            ratio_threshold = float(self.ratio_entry.get())
        except ValueError:
            messagebox.showerror(
                "Input Error",
                "Please fill cell, nucleus, and ratio parameters first."
            )
            return

        if not self.cell_preview_image_path or not self.nucleus_preview_image_path:
            messagebox.showerror(
                "Input Error",
                "Please load both a cellbody example image and a nucleus example image first."
            )
            return

        watershed_mode = self.watershed_mode_var.get()
        ratio_preview_window = tk.Toplevel(self.root)
        ratio_preview = RatioPreviewWindow(
            ratio_preview_window,
            scale,
            intensity_threshold_cell,
            min_size_cell,
            max_size_cell,
            intensity_threshold_nucleus,
            min_size_nucleus,
            max_size_nucleus,
            circularity_threshold=float(self.circularity_threshold_entry.get()),
            watershed_mode=watershed_mode,
            watershed_sigma=float(self.watershed_sigma_entry.get()) if watershed_mode != "No Watershed" else self.watershed_sigma,
            watershed_min_distance=int(self.watershed_min_distance_entry.get()) if watershed_mode != "No Watershed" else self.watershed_min_distance,
            watershed_h_value=float(self.watershed_h_value_entry.get()) if watershed_mode != "No Watershed" else self.watershed_h_value,
            ratio_threshold=ratio_threshold,
            initial_cell_file_path=self.cell_preview_image_path,
            initial_nucleus_file_path=self.nucleus_preview_image_path,
        )
        self.root.wait_window(ratio_preview_window)

        set_entry_value(self.ratio_entry, ratio_preview.ratio_threshold)

        if ratio_preview.cell_file_path:
            self.cell_preview_image_path = ratio_preview.cell_file_path
        if ratio_preview.nucleus_file_path:
            self.nucleus_preview_image_path = ratio_preview.nucleus_file_path

    def run(self):
        self.root.mainloop()
        return self.parameters




# Application control workflow

def collect_image_pairs(cell_folders, nucleus_folders):
    if len(cell_folders) != len(nucleus_folders):
        raise ValueError("The number of cell and nucleus folders must be the same")

    cell_image_files = list_image_files(cell_folders)
    nucleus_image_files = list_image_files(nucleus_folders)

    return list(zip(cell_image_files, nucleus_image_files))


def select_input_folders():
    cell_folders = select_directories_tk("cell")
    nucleus_folders = select_directories_tk("nucleus")
    return cell_folders, nucleus_folders


def main():
    show_loading_screen_and_load()
    parameter_gui = ParameterGUI(DEFAULT_PARAMETERS)
    parameters = parameter_gui.run()

    if parameters is None:
        print("Parameter dialog was closed without saving.")
        return

    cell_folders, nucleus_folders = select_input_folders()

    if not cell_folders or not nucleus_folders:
        print("No folders selected.")
        return

    image_paths_pairs = collect_image_pairs(cell_folders, nucleus_folders)

    save_parameters_to_folders(cell_folders, parameters)
    run_batch_processing(image_paths_pairs, parameters)
    aggregate_csvs(cell_folders)
    aggregate_non_circular_nuclei_csvs(cell_folders)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
