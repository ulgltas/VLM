#!/usr/bin/env python3
# -*- coding: utf8 -*-
# Copyright 2019 Université de Liège
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# VLM run script

import argparse, os
import pythonVLM.VLM_utilities as utils

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs=1, help="Python input file")
    parser.add_argument("--nogui", action="store_true", help="Disable GUI", default=False)
    args = parser.parse_args()


    orig_dir = os.getcwd()                              # Initial working directory

    names = args.file[0].split(os.sep)
    filename = names[-1]                                # Filename of the test case, usually ending on .py
    origin = os.path.split(args.file[0])[0]             # Path from the working directory to the file
    workspace_name = os.path.splitext(args.file[0])[0]  # Workspace name including filename but not .py

    utils.add_path(origin, True)

    workspace = utils.create_workspace(workspace_name, origin)
    print(("Workspace: {}".format(workspace)))

    

    filename = os.path.join(orig_dir, origin, filename) # Absolute path to the file
    if os.path.isfile(filename):
        print(("Executing file {}".format(filename)))
        exec(compile(open(filename, "rb").read(), filename, 'exec'), globals(), locals())
    else:
        raise Exception("File not found: {}".format(filename))


if __name__ == "__main__":
    main()
