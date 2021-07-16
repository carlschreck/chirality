#!/bin/bash

python visualisecells.py $1
/usr/bin/ffmpeg -i myfile%04d.png -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2:color=white" -vcodec libx264 -pix_fmt yuv420p $2
rm myfile*.png
