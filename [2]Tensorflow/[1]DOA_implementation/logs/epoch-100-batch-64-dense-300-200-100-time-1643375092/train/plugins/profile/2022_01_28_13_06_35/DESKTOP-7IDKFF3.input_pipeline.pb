	Gx$h@Gx$h@!Gx$h@	p}C??#*@p}C??#*@!p}C??#*@"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$Gx$h@ŏ1w-!??Ar?????@Y?b?=y??*	     H?@2O
Iterator::Model::BatchV23ı.n???!j:???X@)A??ǘ???1?*?3ڽU@:Preprocessing2Y
!Iterator::Model::BatchV2::Shuffle??z6?>??!?{??6`(@)?z6?>??1?{??6`(@:Preprocessing2F
Iterator::Modelףp=
???!      Y@) ?o_?y?1????v??:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
both?Your program is MODERATELY input-bound because 13.1% of the total step time sampled is waiting for input. Therefore, you would need to reduce both the input time and other time.no*moderate2s3.1 % of the total step time sampled is spent on 'All Others' time. This could be due to Python execution overhead.9p}C??#*@>Look at Section 3 for the breakdown of input time on the host.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	ŏ1w-!??ŏ1w-!??!ŏ1w-!??      ??!       "      ??!       *      ??!       2	r?????@r?????@!r?????@:      ??!       B      ??!       J	?b?=y???b?=y??!?b?=y??R      ??!       Z	?b?=y???b?=y??!?b?=y??JCPU_ONLYYp}C??#*@b 