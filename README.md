Workflow overview:
1. Make config for illumiprocessor with sample names, adapter seqs for raw reads. conf CANNOT be an RTF file, need to download blank .conf file used in previous illumiprocessor script and overwrite. Some sample name formats from GSL are incompatible with illumiprocessor due to underscores and may need to be manually renamed to fit the following convention: <site><year><two digit sample number>
csv headers must be ordered as follows: plate	well	i7	i5	oligoname	libraryname
```
make_illumiprocessor_config.py
```


