[![ascl:2108.006](https://img.shields.io/badge/ascl-2108.006-blue.svg?colorB=262255)](https://ascl.net/2108.006)
[![DOI](https://img.shields.io/badge/DOI-10.1051%2F00046361%2F202553919-blue)](https://doi.org/10.1051/0004-6361/202553919)

# viper - Velocity and IP EstimatoR

For further information, manual, support visit https://mzechmeister.github.io/viper_RV_pipeline.

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

If you publish results with viper, please acknowledge it by citing [Köhler et al. (2025, A&A, 698, 44)](https://ui.adsabs.harvard.edu/abs/2025A&A...698A..44K) [![PDF](https://upload.wikimedia.org/wikipedia/commons/9/94/PDF_icon_-_grey-red_-_16px.svg)](https://www.aanda.org/articles/aa/pdf/2025/06/aa53919-25.pdf). Lower case and monospace font is preferred, i.e. in LaTeX `{\tt viper}`.
