#!/usr/bin/env bash
# * Copy data from the NAS to CHPC for processing

# * Environment
# ** Remote NAS (source)
r1id="vkoppelmans"
r1ip="eindhoven.synology.me"
r1dir="/volume1/Backup/Thalia/Kladblok/20190830_Koppelmans\
/20240816_NMD_fMRIprep_Tethys/data/20240715_BIDS"

# ** CHPC (target)
r2id="u6012627"
r2ip="notchpeak.chpc.utah.edu"
r2dir="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Kladblok\
/20251114_NMD_fMRI"


# * Copy data from NAS to local
# NAS does not allow scp protocol, and rsync does not allow remote-to-remote
# copying, so I will need to copy the data over from NAS to local, and then from
# local to CHPC.

# ** Copy NAS to local
tdir=$(mktemp -d)
rsync -avmhe ssh ${r1id}@${r1ip}:${r1dir} "${tdir}"

# ** Copy local to CHPC
# * Create output folder
ssh "${r2id}@${r2ip}" mkdir -p "${r2dir}"

# ** Copy over the data
rsync -avmhe ssh "${tdir}/20240715_BIDS" ${r2id}@${r2ip}:${r2dir}
