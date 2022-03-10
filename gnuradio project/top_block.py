#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Top Block
# Generated: Thu Mar 10 22:55:58 2022
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import forms
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import time
import wx


class top_block(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Top Block")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 20e6
        self.power = power = 1
        self.center_fre = center_fre = 2.4e9
        self.bandwidth = bandwidth = 20e6

        ##################################################
        # Blocks
        ##################################################
        _power_sizer = wx.BoxSizer(wx.VERTICAL)
        self._power_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_power_sizer,
        	value=self.power,
        	callback=self.set_power,
        	label='power',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._power_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_power_sizer,
        	value=self.power,
        	callback=self.set_power,
        	minimum=0,
        	maximum=100,
        	num_steps=100,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_power_sizer)
        self.wxgui_scopesink2_0 = scopesink2.scope_sink_c(
        	self.GetWin(),
        	title='Scope Plot',
        	sample_rate=samp_rate,
        	v_scale=0,
        	v_offset=0,
        	t_scale=0,
        	ac_couple=False,
        	xy_mode=False,
        	num_inputs=1,
        	trig_mode=wxgui.TRIG_MODE_AUTO,
        	y_axis_label='Counts',
        )
        self.Add(self.wxgui_scopesink2_0.win)
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_c(
        	self.GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title='FFT Plot',
        	peak_hold=False,
        )
        self.Add(self.wxgui_fftsink2_0.win)
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("addr0=192.168.1.23", "")),
        	uhd.stream_args(
        		cpu_format="sc16",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0.set_clock_rate(200e6, uhd.ALL_MBOARDS)
        self.uhd_usrp_source_0.set_clock_source('external', 0)
        self.uhd_usrp_source_0.set_time_source('external', 0)
        self.uhd_usrp_source_0.set_subdev_spec('A:0', 0)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_0.set_center_freq(center_fre, 0)
        self.uhd_usrp_source_0.set_gain(power, 0)
        self.uhd_usrp_source_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_0.set_bandwidth(bandwidth, 0)
        self.blocks_interleaved_short_to_complex_0 = blocks.interleaved_short_to_complex(True, False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_int*1, '/home/lab2/data/wifi_data1.bin', False)
        self.blocks_file_sink_0.set_unbuffered(False)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_interleaved_short_to_complex_0, 0), (self.wxgui_fftsink2_0, 0))
        self.connect((self.blocks_interleaved_short_to_complex_0, 0), (self.wxgui_scopesink2_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.blocks_interleaved_short_to_complex_0, 0))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)
        self.wxgui_fftsink2_0.set_sample_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)

    def get_power(self):
        return self.power

    def set_power(self, power):
        self.power = power
        self._power_slider.set_value(self.power)
        self._power_text_box.set_value(self.power)
        self.uhd_usrp_source_0.set_gain(self.power, 0)

        self.uhd_usrp_source_0.set_gain(self.power, 1)

        self.uhd_usrp_source_0.set_gain(self.power, 2)

        self.uhd_usrp_source_0.set_gain(self.power, 3)

        self.uhd_usrp_source_0.set_gain(self.power, 4)

        self.uhd_usrp_source_0.set_gain(self.power, 5)

        self.uhd_usrp_source_0.set_gain(self.power, 6)

        self.uhd_usrp_source_0.set_gain(self.power, 7)


    def get_center_fre(self):
        return self.center_fre

    def set_center_fre(self, center_fre):
        self.center_fre = center_fre
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 0)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 1)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 2)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 3)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 4)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 5)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 6)
        self.uhd_usrp_source_0.set_center_freq(self.center_fre, 7)

    def get_bandwidth(self):
        return self.bandwidth

    def set_bandwidth(self, bandwidth):
        self.bandwidth = bandwidth
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 0)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 1)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 2)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 3)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 4)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 5)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 6)
        self.uhd_usrp_source_0.set_bandwidth(self.bandwidth, 7)


def main(top_block_cls=top_block, options=None):

    tb = top_block_cls()
    tb.Start(True)
    tb.Wait()


if __name__ == '__main__':
    main()
