<a href="https://ascl.net/2108.006"><img src="https://img.shields.io/badge/ascl-2108.006-blue.svg?colorB=262255" alt="ascl:2108.006" /></a>

# viper - Velocity and IP EstimatoR

```
git clone https://github.com/mzechmeister/viper.git
```

Create shortcuts:
```bash
ln -s $PWD/viper/viper.py ~/bin/viper
ln -s $PWD/viper/vpr.py ~/bin/vpr
ln -s $PWD/viper/GUI_viper.py ~/bin/viper_gui
ln -s $PWD/viper/GUI_vpr.py ~/bin/vpr_gui
```

To run:
```
viper "data/TLS/hd189733/*" data/TLS/Deconv/HARPS*fits -oset 19:21 -nset :4
```
This runs from order 19 (inclusive) to 21 (exclusive) for the first 4 observational files.

To analyse the RVs afterwards use:
```
vpr <tag>
```
`<tag>` defaults to `tmp` in `viper` and `vpr`. See `viper -?` for more options.

If you publish results with viper, please acknowledge it by citing its bibcode from https://ui.adsabs.harvard.edu/abs/2021ascl.soft08006Z.
Lower case and monospace font is preferred, i.e. in LaTeX `{\tt viper}`.
