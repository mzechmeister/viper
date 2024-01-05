<a href="https://ascl.net/2108.006"><img src="https://img.shields.io/badge/ascl-2108.006-blue.svg?colorB=262255" alt="ascl:2108.006" /></a>

# viper - Velocity and IP EstimatoR

For further information, manual, support visit https://mzechmeister.github.io/viper_RV_pipeline .

Download viper:
```
git clone https://github.com/mzechmeister/viper.git
```

Download demo data:
```
git clone https://github.com/mzechmeister/viper_demo_data.git
```

Installation (run from viper directory):
```bash
pip install -e .
```

To run:
```
viper "data/TLS/HD189733/*" data/TLS/HD189733_tpl/HARPS*fits -oset 19:21 -nset :4
```
This runs from order 19 (inclusive) to 21 (exclusive) for the first 4 observational files.

To analyse the RVs afterwards use:
```
vpr <tag>
```
`<tag>` defaults to `tmp` in `viper` and `vpr`. See `viper -?` for more options.

If you publish results with viper, please acknowledge it by citing its bibcode from https://ui.adsabs.harvard.edu/abs/2021ascl.soft08006Z.
Lower case and monospace font is preferred, i.e. in LaTeX `{\tt viper}`.
