# Raw Data Readme
## Syncing
### Which directories to sync
In the file `include_list.txt` include start and end of the rsync command executed
by the Makefile command. Directory names should be terminated in /.
### Which files to include.
This is dictated by the `` file. If you want to exclude everything execept what is 
included in the file the file must end in:
\- */* 
to avoid ignoring directories. Since the
exclusion is based on paths and not files the included files should include the path. To include 
traj_comp.xtc files you should add:
\+ *traj_comp.xtc
