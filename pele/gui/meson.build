
MainWindow_source = custom_target(
 'MainWindow.py',
 input : 'MainWindow.ui',
 output : 'MainWindow.py',
 command : [pyuic, '@INPUT@', '-o', '@OUTPUT@'],
 install : true,
 install_dir : pydir/'pele/gui',
)


connect_explorer_ui_source = custom_target(
 'connect_explorer_ui.py',
 input : 'connect_explorer_ui.ui',
 output : 'connect_explorer_ui.py',
 command : [pyuic, '@INPUT@', '-o', '@OUTPUT@'],
 install : true,
 install_dir : pydir/'pele/gui',
)


ui_params_source = custom_target(
 'ui_params.py',
 input : 'ui_params.ui',
 output : 'ui_params.py',
 command : [pyuic, '@INPUT@', '-o', '@OUTPUT@'],
 install : true,
 install_dir : pydir/'pele/gui',
)

subdir('ui')
