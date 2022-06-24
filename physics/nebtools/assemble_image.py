#! /usr/bin/env python
"""Assemble multiple images into one large image"""

from PIL import Image


def assemble_image(image_names, out_name="allin1.png", margin_ratio=None,
                   grid_size=None, row_major=True):
    """
    Assemble images into one large image.

    :param List[str] image_names: file names of images
    :param str out_name: file name of output image
    :param tuple margin_ratio: margin ratios along x and y directions
    :param tuple grid_size: number of rows and columns
    :param bool row_major: whether the grid is filled in row-major order
    :return: None
    :raises ValueError: if grid_size is not large enough to hold the images
    """
    if margin_ratio is None:
        margin_ratio = (0.05, 0.1)
    if grid_size is None:
        grid_size = (1, len(image_names))
    if len(image_names) > (grid_size[0] * grid_size[1]):
        raise ValueError("Grid size too small")

    # Evaluate image widths and heights
    image_ref = Image.open(image_names[0])
    w0, h0 = image_ref.size
    width_scale = w0 + 2 * int(w0 * margin_ratio[0])
    height_scale = h0 + 2 * int(h0 * margin_ratio[1])
    width_total = width_scale * grid_size[1]
    height_total = height_scale * grid_size[0]

    # Assemble images
    image_total = Image.new("RGB", (width_total, height_total), (255, 255, 255))
    margin_x = int(w0 * margin_ratio[0])
    margin_y = int(h0 * margin_ratio[1])
    if row_major:
        for i in range(grid_size[0]):
            dy = margin_y + height_scale * i
            for j in range(grid_size[1]):
                dx = margin_x + width_scale * j
                idx = grid_size[1] * i + j
                image_idx = Image.open(image_names[idx])
                image_total.paste(image_idx, (dx, dy))
    else:
        for j in range(grid_size[1]):
            dx = margin_x + width_scale * j
            for i in range(grid_size[0]):
                dy = margin_y + height_scale * i
                idx = grid_size[0] * j + i
                image_idx = Image.open(image_names[idx])
                image_total.paste(image_idx, (dx, dy))
    image_total.save(out_name)


def main():
    image_names = [f"{i}.png" for i in range(8)]
    assemble_image(image_names, margin_ratio=(0.05, 0.2))


if __name__ == "__main__":
    main()
