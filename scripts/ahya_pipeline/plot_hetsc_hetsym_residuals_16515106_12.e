Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/compare_model.py", line 128, in <module>
    main(args.data, args.model, args.mask, args.fold, args.opt, PTS, args.resid_range)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/compare_model.py", line 103, in main
    model_func = Numerics.make_extrap_log_func(SETTINGS.get_settings(model))
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/SETTINGS.py", line 267, in get_settings
    raise ValueError(
ValueError: Model nickname undefined please check in SETTINGS.py you are using the correct model nickname!
