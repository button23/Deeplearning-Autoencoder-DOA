?	,e?X@,e?X@!,e?X@	?v+??c???v+??c??!?v+??c??"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$,e?X@?"??~j??Ash??|?@Y?sF????*	?????Yk@2v
?Iterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate^?I+??!??e-U?T@)??ͪ????1e?HbT@:Preprocessing2l
5Iterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat46<???!??=TD? @)???Q???1f??l@:Preprocessing2F
Iterator::Model-C??6??!?Z?\~f@)?<,Ԛ?}?1?O?!??
@:Preprocessing2U
Iterator::Model::ParallelMapV2?I+?v?!f??@)?I+?v?1f??@:Preprocessing2Z
#Iterator::Model::ParallelMapV2::Zip??|?5^??!SZ5??W@)??ZӼ?t?1Z?ױ??@:Preprocessing2x
AIterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat::FromTensorF%u?k?!???O?!??)F%u?k?1???O?!??:Preprocessing2?
OIterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice??_?Le?!?iJ?F??)??_?Le?1?iJ?F??:Preprocessing2f
/Iterator::Model::ParallelMapV2::Zip[0]::FlatMap???o_??!s?*"?T@)-C??6Z?1?Z?\~f??:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
device?Your program is NOT input-bound because only 0.6% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no9?v+??c??#You may skip the rest of this page.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	?"??~j???"??~j??!?"??~j??      ??!       "      ??!       *      ??!       2	sh??|?@sh??|?@!sh??|?@:      ??!       B      ??!       J	?sF?????sF????!?sF????R      ??!       Z	?sF?????sF????!?sF????JCPU_ONLYY?v+??c??b Y      Y@q???3V:@"?
device?Your program is NOT input-bound because only 0.6% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.b
`input_pipeline_analyzer (especially Section 3 for the breakdown of input operations on the Host)m
ktrace_viewer (look at the activities on the timeline of each Host Thread near the bottom of the trace view)"T
Rtensorflow_stats (identify the time-consuming operations executed on the CPU_ONLY)"Z
Xtrace_viewer (look at the activities on the timeline of each CPU_ONLY in the trace view)*?
?<a href="https://www.tensorflow.org/guide/data_performance_analysis" target="_blank">Analyze tf.data performance with the TF Profiler</a>*y
w<a href="https://www.tensorflow.org/guide/data_performance" target="_blank">Better performance with the tf.data API</a>2I
=type.googleapis.com/tensorflow.profiler.GenericRecommendation
nono:
Refer to the TF2 Profiler FAQb?26.3367% of Op time on the host used eager execution. Performance could be improved with <a href="https://www.tensorflow.org/guide/function" target="_blank">tf.function.</a>2"CPU: B 