	?O0@?O0@!?O0@	$?P??@$?P??@!$?P??@"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$?O0@?f??j+??A?y?):/@Y?d?`TR??*	???????@2F
Iterator::Model/n????!?_??3FJ@)F????x??1?<)???I@:Preprocessing2v
?Iterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate?/?'??!	Ԣn?j9@)[????<??13?+(?8@:Preprocessing2f
/Iterator::Model::ParallelMapV2::Zip[0]::FlatMap<?R?!???!?Y_]F@)?J?4??1F+E?O3@:Preprocessing2l
5Iterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat46<???!>?#?Qf @)? ?	???1'?h@???:Preprocessing2U
Iterator::Model::ParallelMapV2"??u????!ʈb????)"??u????1ʈb????:Preprocessing2?
OIterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice ?o_?y?!?:??J???) ?o_?y?1?:??J???:Preprocessing2Z
#Iterator::Model::ParallelMapV2::Zipxz?,C??!?Lp̹G@)?~j?t?x?1??8????:Preprocessing2x
AIterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat::FromTensor?????g?!Uz?u???)?????g?1Uz?u???:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
device?Your program is NOT input-bound because only 3.0% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no9$?P??@#You may skip the rest of this page.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	?f??j+???f??j+??!?f??j+??      ??!       "      ??!       *      ??!       2	?y?):/@?y?):/@!?y?):/@:      ??!       B      ??!       J	?d?`TR???d?`TR??!?d?`TR??R      ??!       Z	?d?`TR???d?`TR??!?d?`TR??JCPU_ONLYY$?P??@b 