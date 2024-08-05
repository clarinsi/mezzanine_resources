# 2024-08-02T11:53:05

## Inference:

With 0_testing_inference it was found out that the model runs reliably for audio durations all the way to 200s, after which testing was stopped. This should be repeated on GPU.

## Training: (2024-08-05T09:23:37)

Training now works. Key points:
* in `feature_extractor`, `max_length` parameter should be the lenght of the desired output (50 features per second of duration.)
* We can specify the number of classes in model config (i.e. `Model(num_labels=3)`), but the labels must be one-hot encoded (`[1,0,0]`)
* Padding: not to be called from `feature_extractor(audio_array, ...)`. The outputs of the feature extraction have to be input anew to `feature_extractor.pad(inputs, pad="max_length")`, else it doesn't work.