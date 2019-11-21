#!/usr/bin/env python
# -*- coding: utf8 -*-
# VLM run script

import argparse, os
import pythonVLM.VLM_utilities as utils

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs=1, help="Python input file")
    parser.add_argument("--nogui", action="store_true", help="Disable GUI", default=False)
    args = parser.parse_args()


    vlm_dir = os.path.abspath(os.path.split(__file__)[0])

    names = args.file[0].split(os.sep)
    filename = names[-1]
    name = names[-2]
    origin = os.path.split(args.file[0])[0]

    utils.add_path(os.path.join("VLM",origin))

    workspace = utils.create_workspace(name, origin)
    print("Workspace: {}".format(workspace))

    

    filename = os.path.join(vlm_dir, origin, filename)
    if os.path.isfile(filename):
        print("Executing file {}".format(filename))
        execfile(filename, globals(), locals())
    else:
        raise Exception("File not found: {}".format(filename))


if __name__ == "__main__":
    main()