"""Install the optional pinned RapidNJ source backend explicitly."""

import argparse
from chronoclade.cgmlst.scale_setup import install_rapidnj

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    print(install_rapidnj(args.target, jobs=args.jobs))
