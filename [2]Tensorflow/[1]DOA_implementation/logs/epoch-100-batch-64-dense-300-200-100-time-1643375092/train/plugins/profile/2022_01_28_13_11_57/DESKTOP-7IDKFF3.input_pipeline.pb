	h"lxz%)@h"lxz%)@!h"lxz%)@	?-˥#@?-˥#@!?-˥#@"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$h"lxz%)@6<?R???AD?l???&@Y?d?`TR??*	    ԍ@2[
$Iterator::Model::FiniteTake::BatchV2 c?ZB>??!?X????X@)?sF????1?K?V@:Preprocessing2e
-Iterator::Model::FiniteTake::BatchV2::Shuffle??|?5^???!?????%@)?|?5^???1?????%@:Preprocessing2F
Iterator::Model%??C???!      Y@)U???N@s?1ÜӾ????:Preprocessing2R
Iterator::Model::FiniteTake??7??d??!c,Ar|?X@)U???N@s?1ÜӾ????:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
both?Your program is MODERATELY input-bound because 7.8% of the total step time sampled is waiting for input. Therefore, you would need to reduce both the input time and other time.no*no9?-˥#@>Look at Section 3 for the breakdown of input time on the host.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	6<?R???6<?R???!6<?R???      ??!       "      ??!       *      ??!       2	D?l???&@D?l???&@!D?l???&@:      ??!       B      ??!       J	?d?`TR???d?`TR??!?d?`TR??R      ??!       Z	?d?`TR???d?`TR??!?d?`TR??JCPU_ONLYY?-˥#@b 