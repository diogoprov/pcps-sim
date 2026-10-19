# Raw data (immutable)

- `portik2023_timetree.tre`: species-level amphibian timetree of Portik, Streicher &
  Wiens 2023, *Mol. Phylogenet. Evol.*, doi:10.1016/j.ympev.2023.107907 (5,326 tips,
  including non-anuran outgroups). Copied on 2026-10-17 from the author's compendium
  `~/GitHub/brazilian-tadpole-knowledge-gaps/data/portik2023_timetree.tre`.
  SHA-256 d871f3643575709e846f167ce297834a63afe21743296902d1548db4b4691e90.
  Single (point-estimate) tree. It is byte-identical to
  `TreePL-Rooted_Anura_bestTree.tre` in the authors' OSF component "Final Trees"
  (https://osf.io/qgcf2/, part of https://osf.io/3cbsh/; licence GPL-3.0,
  copyright Daniel Portik 2023; checked 2026-10-18). The same component holds
  the time-calibrated bootstrap trees (`TreePL-Rooted_Anura_bootstraps.tre`) and
  the node confidence intervals; there is no Bayesian posterior.

  The tree file is not redistributed in this repository (it is listed in
  `.gitignore`). Download it from https://osf.io/download/akupm/, save it as
  `data/raw/portik2023_timetree.tre` and check the SHA-256 above before running
  the pipeline.

- `portik2023_timetree_bootstraps.tre`: the 100 time-calibrated bootstrap trees
  (`TreePL-Rooted_Anura_bootstraps.tre`, OSF file id 64ea0805209bae05d5e93844,
  same component and licence as above). Downloaded on 2026-10-19 from
  https://files.osf.io/v1/resources/qgcf2/providers/osfstorage/64ea0805209bae05d5e93844?action=download.
  SHA-256 b608766078c2570cf76bd342d0c4232c68cb3de57199c9270245289b891212db.
  Not redistributed here (listed in `.gitignore`).
