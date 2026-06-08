#!/usr/bin/env python3
import sys
import csv
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract trait metadata from fastGWA CSV based on GWA_SUMSTA"
    )
    parser.add_argument(
        "--gwa-sumsta",
        required=True,
        help="GWA summary statistics file name (e.g. xxx.v1.1.fastGWA.gz)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gwa = args.gwa_sumsta

    if gwa.endswith(".v1.1.fastGWA.gz"):
        csv_file = "inputs/fastgwa_UKB_impute_v1.1.csv"
        mode = "v1.1"
    elif gwa.endswith(".v1.0.fastGWA.gz"):
        csv_file = "inputs/fastgwa_UKB_binary_v1.11.csv"
        mode = "v1.0"
    else:
        sys.exit(f"Unrecognized GWA_SUMSTA format: {gwa}")

    found = False

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            url = row["URL"]
            if os.path.basename(url) != gwa:
                continue

            found = True

            if mode == "v1.1":
                data_type_raw = row["Data_type"]
                N = row["N"]

                if data_type_raw in ("Continuous", "Ordered_Categorical"):
                    DATA_TYPE = "quant"
                    N_CASE = "NA"
                    N_CONTROL = "NA"
                elif data_type_raw == "Binary":
                    DATA_TYPE = "cc"
                    N_CASE = row["Ncase"]
                    N_CONTROL = str(int(row["N"]) - int(row["Ncase"]))
                else:
                    sys.exit(f"Unknown Data_type: {data_type_raw}")

            else:  # v1.0
                DATA_TYPE = "cc"
                N = row["N"]
                N_CASE = row["N_case"]
                N_CONTROL = row["N_control"]

            # bash で eval できる形式
            print(f'TRAIT_TYPE="{DATA_TYPE}"')
            print(f'N="{N}"')
            print(f'N_CASE="{N_CASE}"')
            print(f'N_CONTROL="{N_CONTROL}"')
            break

    if not found:
        sys.exit(f"{gwa} not found in {csv_file}")


if __name__ == "__main__":
    main()
