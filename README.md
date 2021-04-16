# VIPER - Velocity and IP EstimatoR

```
git clone https://github.com/mzechmeister/VIPER.git
```

Create shortcuts:
```bash
ln -s $PWD/viper.py ~/bin/viper
ln -s $PWD/vpr.py ~/bin/vpr
```

To run:
```
viper "data/TLS/hd189733/*" data/TLS/Deconv/HARPS*fits -oset 19:21 -nset :4
```
This runs from order 19(inclusive) to 21(exclusive) for the first 4 observational files.

To analyse the RVs afterwards use:
```
vpr <tag>
```
`<tag>` defaults to `tmp` in `viper` and `vpr`. See `viper -?` for more options.
