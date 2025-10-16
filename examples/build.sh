#!/bin/bash

meson setup --wipe build  && meson compile -vvv -C build

for so in $(find build/new_potential -iname "*.so"); do
  cp ${so} ${so#*/};
done

cp build/gui/* gui/

