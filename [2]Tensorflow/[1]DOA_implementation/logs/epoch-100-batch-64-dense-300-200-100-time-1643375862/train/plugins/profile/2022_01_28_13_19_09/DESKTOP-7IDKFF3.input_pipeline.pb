	?[ Aa5@?[ Aa5@!?[ Aa5@	????????????!??????"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$?[ Aa5@o??ʡ??A?&SC5@Y?J?4??*	23333?d@2v
?Iterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate?46<??!?????eQ@))\???(??1U)?P@:Preprocessing2l
5Iterator::Model::ParallelMapV2::Zip[1]::ForeverRepeatݵ?|г??!?^y?R?.@)'???????135???*@:Preprocessing2Z
#Iterator::Model::ParallelMapV2::Zip?A`??"??!_y?R??V@)	?^)ˀ?1?-|D??@:Preprocessing2F
Iterator::Model?<,Ԛ???!
5?kE?!@)????Mb??1?K???@:Preprocessing2U
Iterator::Model::ParallelMapV2F%u?{?!l??v@)F%u?{?1l??v@:Preprocessing2?
OIterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice?J?4q?!r??y@)?J?4q?1r??y@:Preprocessing2x
AIterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat::FromTensor???_vOn?!???	@)???_vOn?1???	@:Preprocessing2f
/Iterator::Model::ParallelMapV2::Zip[0]::FlatMap:??H???!??o?޳Q@)????Mb`?1?K?????:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
device?Your program is NOT input-bound because only 0.2% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no9??????#You may skip the rest of this page.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	o??ʡ??o??ʡ??!o??ʡ??      ??!       "      ??!       *      ??!       2	?&SC5@?&SC5@!?&SC5@:      ??!       B      ??!       J	?J?4???J?4??!?J?4??R      ??!       Z	?J?4???J?4??!?J?4??JCPU_ONLYY??????b 