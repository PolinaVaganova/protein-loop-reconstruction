import os
import re


def prep_mtz(mtz_in, mtz_out):
    """
    Prepare an MTZ file by doubling the unit cell dimensions and indices.

    Args:
        mtz_in (str): Input MTZ file path.
        mtz_out (str): Output MTZ file path.

    Raises:
        FileExistsError: If the output file already exists.
        AssertionError: If the MTZ file spacegroup is not 'P 1'.
        RuntimeError: If the MTZ file is missing the R-FREE-FLAG column.

    Returns:
        None
    """
    # Avoid overwriting the existing file
    if os.path.isfile(mtz_out):
        raise FileExistsError(f"Output file '{mtz_out}' already exists.")

    # Read the input MTZ file
    mtz = gemmi.read_mtz_file(mtz_in)

    # Ensure the MTZ file has the required spacegroup 'P 1'
    assert mtz.spacegroup.hm == "P 1", "MTZ file spacegroup must be 'P 1'."

    # Regular expression pattern to search for the R-FREE-FLAG column
    R_FREE_FLAG = None
    r_free_pattern_string = r"(r.*free)|(free.*r)"
    r_free_pattern = re.compile(r_free_pattern_string, flags=re.IGNORECASE)

    # Find the R-FREE-FLAG column in the MTZ file
    for column in mtz.columns:
        if r_free_pattern.search(column.label):
            R_FREE_FLAG = column
            break

    # Raise an error if the R-FREE-FLAG column is not found
    if R_FREE_FLAG is None:
        raise RuntimeError(
            f"MTZ file missing R-FREE-FLAG column "
            f"(pattern: `{r_free_pattern_string}`)"
            f"\nPresent columns: {[column.label for column in mtz.columns]}"
        )

    # Get the required columns from the MTZ file
    H, K, L, FOBS, SIGMA_FOBS = [
        mtz.column_with_label(label) for label in ("H", "K", "L", "FOBS", "SIGFOBS")
    ]

    # Create a default SIGMA_FOBS column if it is missing
    if SIGMA_FOBS is None:
        SIGMA_FOBS = np.ones(len(FOBS))

    # Create a list to store the doubled data
    data = []

    # Double the indices and add data to the list
    for h, k, l, fobs, sigma, r_flag in zip(H, K, L, FOBS, SIGMA_FOBS, R_FREE_FLAG):
        data.append([h * 2, k * 2, l * 2, fobs, sigma, r_flag])

    # Convert the data list to a numpy array
    data = np.array(data)

    # Create a new MTZ object with doubled unit cell dimensions
    mtzx2 = gemmi.Mtz()
    mtzx2.spacegroup = mtz.spacegroup
    double_uc = mtz.get_cell()
    mtzx2.set_cell_for_all(
        gemmi.UnitCell(
            double_uc.a * 2,
            double_uc.b * 2,
            double_uc.c * 2,
            double_uc.alpha,
            double_uc.beta,
            double_uc.gamma,
        )
    )

    # Add a new dataset and columns to the MTZ object
    mtzx2.add_dataset("Doubled Indices")
    mtzx2.add_column("H", "H")
    mtzx2.add_column("K", "H")
    mtzx2.add_column("L", "H")
    mtzx2.add_column("F", "F")
    mtzx2.add_column("SIGF", "Q")
    mtzx2.add_column("R-FREE", "I")

    # Set the data in the MTZ object
    mtzx2.set_data(data)

    # Write the MTZ object to the output file
    mtzx2.write_to_file(mtz_out)


if __name__ == "__main__":
    # Example usage
    # prep_mtz("input.mtz", "output.mtz")
    pass
