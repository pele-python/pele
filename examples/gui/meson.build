
_blj_dialog_source = custom_target(
 '_blj_dialog.py',
 input : '_blj_dialog.ui',
 output : '_blj_dialog.py',
 command : [pyuic, '@INPUT@', '-o', '@OUTPUT@'],
 install : true,
 install_dir : pydir/'gui'
)
