#!/usr/bin/env python
# * Libraries
import argparse
from glob import glob
import json
from shutil import copy2


# * Import arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Populate fieldmap json files intended for fields'
    )
    parser.add_argument('--sid', required=True, help='Subject ID.')
    parser.add_argument('--ses', required=True, help='Session ID.')
    parser.add_argument('--path', required=True, help='Root folder containing subject folders')
    # ** Parse arguments
    args = parser.parse_args()


# * Create class for patching
class fixjson:
    def __init__(self, sid, ses, path):

        # ** Store input arguments
        self.subid = str(sid)
        self.sesid = str(ses)
        self.base = path

        # ** Build paths
        self.idir = self.base + "/sub-" + self.subid + "/ses-" + self.sesid
        self.idir_func = self.idir + "/func"
        self.idir_fmap = self.idir + "/fmap"

    # ** Locate funtional data
    def locate_nii(self, sequence):
        # *** Find the nifti data paths
        self.niis = sorted(glob(self.idir_func + f"/*{sequence}*nii.gz"))
        # *** Remove root + subject path from these paths
        self.niis = [
            x.replace(self.base + "/sub-" + self.subid + "/", "")
            for x in self.niis
        ]
        return self.niis

    # ** Locate field map json files
    # This should bring up exactly two files: AP and PA for a particular sequence
    def locate_json(self, sequence):
        # *** Find the json data paths
        self.json = sorted(glob(self.idir_fmap + f"/*{sequence}*json"))
        return self.json

    # ** Loop over json files and inject 'intended for' relative path
    def doinject(self, sequence):

        # *** Loop over the json AP/PA images
        for my_json in self.locate_json(sequence):

            # **** Create backup
            # file_bkp = my_json.replace("json", "json_bkp")
            # copy2(my_json, file_bkp)

            # **** Load json file
            with open(my_json) as json_file:
                self.tmp_json = json.load(json_file)

                # ***** Announce
                print(f"Populating IntendedFor field of {my_json.split('/')[-1]} with {self.locate_nii(sequence)[0]}")

                # ***** Inject fMRI path into 'intended for' field
                self.tmp_json['IntendedFor'] = self.locate_nii(sequence)

            # **** Save the edited json over the original file
            with open(my_json, "w") as write_file:
                json.dump(self.tmp_json, write_file, indent=4)


# * Apply patch
# ** Test
# my_fix_object = fixjson("NMD001", "01", "/datadisk/tmp/NMD")
# for sequence in ['handdom', 'rest']:
#     my_fix_object.doinject(sequence)

# ** Production
my_fix_object = fixjson(args.sid, args.ses, args.path)
for sequence in ['handdom', 'rest']:
    my_fix_object.doinject(sequence)
