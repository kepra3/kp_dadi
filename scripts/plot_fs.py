import argparse
import os
import numpy as np
import dadi
import matplotlib
matplotlib.use('Agg')
import pylab
import copy


def apply_mask(fs, mask):
    """Apply masking based on the specified type."""
    if mask == "low" and len(fs.sample_sizes) == 2:
        indices = [(1, 0), (0, 1), (2, 0), (0, 2), (1, 1)]
        for i, j in indices:
            fs.mask[i, j] = True
    elif mask == "mid" and len(fs.sample_sizes) == 1:
        mid = int(fs.sample_sizes[0] / 2)
        fs.mask[mid] = True
    else:
        print("No masking applied")


def calculate_statistic(fs):
    """Calculate the appropriate statistic based on sample size dimensions."""
    if len(fs.sample_sizes) == 2:
        # maybe have a look at using fs.Zengs_E()
        return "FST", fs.Fst()
    elif len(fs.sample_sizes) == 1:
        return "TajimasD", fs.Tajima_D()
    elif len(fs.sample_sizes) == 3:
        return "None", 0.0
    else:
        raise ValueError("Unsupported dimensionality for the site frequency spectrum")


def write_statistics(stats_out_name, snps, fs, statistic_name, statistic_value):
    """Write statistics to the output file."""
    with open(stats_out_name, 'a') as stats_out:
        if stats_out.tell() == 0:
            print('Creating a new file')
            stats_out.write("Pop\tSample sizes\tSum of SFS\t{}\n".format(statistic_name))
        stats_out.write(
            f"{snps}\t{fs.sample_sizes}\t{np.around(fs.S(), 2)}\t{np.around(statistic_value, 2)}\n"
        )



def project_fs(fs, method):
    """Project the frequency spectrum if requested."""
    extras = ""
    if method == "projection":
        if len(fs.sample_sizes) == 1:
            project_int = int(fs.sample_sizes[0] * 0.8)
            fs = fs.project([project_int])
        else:
            project_int = [int(x * 0.8) for x in fs.sample_sizes]
            fs = fs.project(project_int)
        extras += "project"
    return fs, extras


def fold_fs(fs, fold):
    """Fold the frequency spectrum if requested."""
    extras = ""
    if fold == "folded":
        fs = fs.fold()
        extras += "fold"
    return fs, extras

def plot_spectrum(fs, snps, outpath, extras):
    """Plot the site frequency spectrum."""
    # Plotting
    fig = pylab.figure(figsize=(5, 5))
    if len(fs.sample_sizes) == 1:
        dadi.Plotting.plot_1d_fs(fs)
    elif len(fs.sample_sizes) > 1:
        if len(fs.sample_sizes) == 2:
            colour_map = copy.copy(pylab.cm.get_cmap("hsv"))
            dadi.Plotting.plot_single_2d_sfs(fs, vmin=1, cmap=colour_map)
        elif len(fs.sample_sizes) == 3:
            from mayavi import mlab
            # Mask monomorphic bins (optional)
            fs.mask[0, :, :] = True
            fs.mask[:, 0, :] = True
            fs.mask[:, :, 0] = True

            # Convert masked values to zero or NaN (choose depending on your plotting preference)
            data = fs.data.copy()
            data[fs.mask] = 0  # or np.nan if using mlab.points3d

            # Optionally apply log scaling
            data = np.log10(data + 1)

            # Create a 3D grid
            x, y, z = np.mgrid[0:fs.shape[0], 0:fs.shape[1], 0:fs.shape[2]]

            # Plot using volume rendering
            mlab.figure(size=(800, 600))
            vol = mlab.pipeline.scalar_field(x, y, z, data)
            mlab.pipeline.volume(vol, vmin=data.min(), vmax=data.max())

            mlab.axes(xlabel='Group2', ylabel='Group3', zlabel='Group4')
            mlab.colorbar(title='log10(Count+1)', orientation='vertical')
            mlab.title('3D Joint SFS (3 Populations)')
            mlab.show()  

    # Save processed .fs file
    out_dir = os.path.join("../data/fs", outpath) if outpath else "../data/fs"
    os.makedirs(out_dir, exist_ok=True)
    fs_file_path = os.path.join(out_dir, f"{snps}_{extras}.fs")
    fs.to_file(fs_file_path)
    print(f"Saved fs to {fs_file_path}")

    # Save plot
    fig.tight_layout()
    plot_file_path = os.path.join("../plots/spectra/", f"{snps}_{extras}.png")
    fig.savefig(plot_file_path, dpi=300)
    print(f"Saved plot to {plot_file_path}")
    pylab.close(fig)




def main(filepath, fold, method, mask, outpath):
    """Main script functionality."""
    # Extract SNPs name from the file path
    snps = os.path.splitext(os.path.basename(filepath))[0]

    # Import spectrum
    fs = dadi.Spectrum.from_file(filepath)

    # Apply masking
    apply_mask(fs, mask)
    if mask != "no":
        extras_mask = f"mask{mask}"

    # Project spectrum if requested
    fs, extras_proj = project_fs(fs, method)

    # Fold spectrum if requested
    fs, extras_fold = fold_fs(fs, fold)

    # Combine extras for file naming
    extras = extras_proj + extras_fold + extras_mask

    # Calculate statistics
    statistic_name, statistic = calculate_statistic(fs)

    # Write statistics to file
    stats_out_name = "../results/sfs_stats.txt"
    write_statistics(stats_out_name, snps, fs, statistic_name, statistic)

    # Plot spectrum
    plot_spectrum(fs, snps, mask, outpath, extras)


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(
        prog="Plotting fs",
        description="A script for plotting fs from .fs file.",
        usage="%(prog)s [options] <filepath> <fold> <method> <mask>"
    )

    # Required positional arguments
    parser.add_argument(
        "filepath",
        help="Path to the data .fs file."
    )
    parser.add_argument(
        "fold",
        help="If you want your spectrum to be folded (e.g., 'folded'). "
             "Only use on unfolded spectrum!"
    )
    parser.add_argument(
        "method",
        help="Method to use for projection or subsample (e.g., 'projection, "
             "'subsample' or 'no')"
    )
    parser.add_argument(
        "mask",
        help="Type of masking to use (e.g., 'mid', 'low', or 'no')."
    )

    # Optional arguments
    parser.add_argument(
        "-o", "--outpath",
        help="Optional subdirectory for saving processed .fs file and plots. "
             "Defaults to '../data/fs/'."
    )

    # Parse arguments
    args = parser.parse_args()

    # Call the main function
    main(args.filepath, args.fold, args.method, args.mask, args.outpath)
