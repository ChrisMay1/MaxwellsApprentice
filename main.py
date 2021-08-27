#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Morse Code Transmit
# Author: C May
# Generated: Thu Aug 26 19:56:52 2021
##################################################

from distutils.version import StrictVersion

from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from grc_gnuradio import blks2 as grc_blks2
from optparse import OptionParser
import osmosdr
import sip
import sys
import time

class cw_transmit(gr.top_block):
    def __init__(self):
        gr.top_block.__init__(self, "Morse Code Transmit")
        ##################################################
        # Variables
        ##################################################
        self.signal_source = signal_source = 0
        self.samp_rate = samp_rate = int(44.1e3)
        self.audio_rate = audio_rate = int(44.1e3)
        self.rf_samp_rate = rf_samp_rate = 2.205e6
        self.rf_gain = rf_gain = 10
        self.if_gain = if_gain = 40
        self.freq = freq = 3.55e6
        self.bb_gain = bb_gain = 20

        ##################################################
        # Blocks
        ##################################################

        # HackRF
        self.osmosdr_sink_0 = osmosdr.sink( args="numchan=" + str(1) + " " + '' )
        self.osmosdr_sink_0.set_sample_rate(rf_samp_rate)
        self.osmosdr_sink_0.set_center_freq(freq, 0)
        self.osmosdr_sink_0.set_freq_corr(0, 0)
        self.osmosdr_sink_0.set_gain(rf_gain, 0)
        self.osmosdr_sink_0.set_if_gain(if_gain, 0)
        self.osmosdr_sink_0.set_bb_gain(bb_gain, 0)
        self.osmosdr_sink_0.set_antenna('', 0)
        self.osmosdr_sink_0.set_bandwidth(0, 0)
        
        # Wave File Source - 52MHz
        self.blocks_wavfile_source_52MHz = blocks.wavfile_source('/home/pi/Downloads/52MHz_JM.wav', True)
        self.analog_wfm_tx_52MHz = analog.wfm_tx(
        	audio_rate=samp_rate,
        	quad_rate=samp_rate,
        	tau=75e-6,
        	max_dev=5e3,
        	fh=-1.0,
        )
        self.rational_resampler_52MHz = filter.rational_resampler_ccc(
                interpolation=50,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        
        # Wave File Source - 145MHz
        self.blocks_wavfile_source_145MHz = blocks.wavfile_source('/home/pi/Downloads/145MHz_JM.wav', True)
        self.analog_wfm_tx_145MHz = analog.wfm_tx(
        	audio_rate=samp_rate,
        	quad_rate=samp_rate,
        	tau=75e-6,
        	max_dev=5e3,
        	fh=-1.0,
        )
        self.rational_resampler_145MHz = filter.rational_resampler_ccc(
                interpolation=50,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        
        # Wave File Source  - 225MHz
        self.blocks_wavfile_source_225MHz = blocks.wavfile_source('/home/pi/Downloads/225MHz_JM.wav', True)
        self.analog_wfm_tx_225MHz = analog.wfm_tx(
        	audio_rate=samp_rate,
        	quad_rate=samp_rate,
        	tau=75e-6,
        	max_dev=5e3,
        	fh=-1.0,
        )
        self.rational_resampler_225MHz = filter.rational_resampler_ccc(
                interpolation=50,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        
        # Wave File Source  - 435MHz
        self.blocks_wavfile_source_435MHz = blocks.wavfile_source('/home/pi/Downloads/435MHz_JM.wav', True)
        self.analog_wfm_tx_435MHz = analog.wfm_tx(
        	audio_rate=samp_rate,
        	quad_rate=samp_rate,
        	tau=75e-6,
        	max_dev=5e3,
        	fh=-1.0,
        )
        self.rational_resampler_435MHz = filter.rational_resampler_ccc(
                interpolation=50,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        
        # Wave File Source  - 905MHz
        self.blocks_wavfile_source_905MHz = blocks.wavfile_source('/home/pi/Downloads/905MHz_JM.wav', True)
        self.analog_wfm_tx_905MHz = analog.wfm_tx(
        	audio_rate=samp_rate,
        	quad_rate=samp_rate,
        	tau=75e-6,
        	max_dev=5e3,
        	fh=-1.0,
        )
        self.rational_resampler_905MHz = filter.rational_resampler_ccc(
                interpolation=50,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        
        # Morse Code Source
        self.blocks_vector_source = blocks.vector_source_c((0,0,0,1,1,0,1,0,1,1,0,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,), True, 1, [])
        self.blocks_repeat = blocks.repeat(gr.sizeof_gr_complex*1, int(2.56e3))
        self.analog_sig_source = analog.sig_source_c(samp_rate, analog.GR_COS_WAVE, 1000, 1, 0)
        self.blocks_multiply = blocks.multiply_vcc(1)
        self.rational_resampler_morse = filter.rational_resampler_ccc(
                interpolation=50,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        
        
        # Selector Block
        self.blks2_selector = grc_blks2.selector(
        	item_size=gr.sizeof_gr_complex*1,
        	num_inputs=6,
        	num_outputs=1,
        	input_index = signal_source,
        	output_index=0,
        )
        
        

        ##################################################
        # Connections
        ##################################################
        
        
        # Morse Code Connections
        self.connect((self.blocks_vector_source, 0), (self.blocks_repeat, 0))
        self.connect((self.blocks_repeat, 0), (self.blocks_multiply, 0))
        self.connect((self.analog_sig_source, 0), (self.blocks_multiply, 1))
        self.connect((self.blocks_multiply, 0), (self.rational_resampler_morse, 0))
        self.connect((self.rational_resampler_morse, 0), (self.blks2_selector, 0))
        
        
        # FM Transmit Connections
        self.connect((self.blocks_wavfile_source_52MHz, 0), (self.analog_wfm_tx_52MHz, 0))
        self.connect((self.analog_wfm_tx_52MHz, 0), (self.rational_resampler_52MHz, 0))
        self.connect((self.rational_resampler_52MHz, 0), (self.blks2_selector, 1))
        
        self.connect((self.blocks_wavfile_source_145MHz, 0), (self.analog_wfm_tx_145MHz, 0))
        self.connect((self.analog_wfm_tx_145MHz, 0), (self.rational_resampler_145MHz, 0))
        self.connect((self.rational_resampler_145MHz, 0), (self.blks2_selector, 2))
        
        self.connect((self.blocks_wavfile_source_225MHz, 0), (self.analog_wfm_tx_225MHz, 0))
        self.connect((self.analog_wfm_tx_225MHz, 0), (self.rational_resampler_225MHz, 0))
        self.connect((self.rational_resampler_225MHz, 0), (self.blks2_selector, 3))
        
        self.connect((self.blocks_wavfile_source_435MHz, 0), (self.analog_wfm_tx_435MHz, 0))
        self.connect((self.analog_wfm_tx_435MHz, 0), (self.rational_resampler_435MHz, 0))
        self.connect((self.rational_resampler_435MHz, 0), (self.blks2_selector, 4))
        
        self.connect((self.blocks_wavfile_source_905MHz, 0), (self.analog_wfm_tx_905MHz, 0))
        self.connect((self.analog_wfm_tx_905MHz, 0), (self.rational_resampler_905MHz, 0))
        self.connect((self.rational_resampler_905MHz, 0), (self.blks2_selector, 5))
        
        
        # Connect to Hardware
        self.connect((self.blks2_selector, 0), (self.osmosdr_sink_0, 0))

    def get_signal_source(self):
        return self.signal_source

    def set_signal_source(self, signal_source):
        self.signal_source = signal_source
        self.blks2_selector.set_input_index(int(self.signal_source))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_sampling_freq(self.samp_rate)

    def get_rf_samp_rate(self):
        return self.rf_samp_rate

    def set_rf_samp_rate(self, rf_samp_rate):
        self.rf_samp_rate = rf_samp_rate
        self.osmosdr_sink_0.set_sample_rate(self.rf_samp_rate)

    def get_rf_gain(self):
        return self.rf_gain

    def set_rf_gain(self, rf_gain):
        self.rf_gain = rf_gain
        self.osmosdr_sink_0.set_gain(self.rf_gain, 0)

    def get_if_gain(self):
        return self.if_gain

    def set_if_gain(self, if_gain):
        self.if_gain = if_gain
        self.osmosdr_sink_0.set_if_gain(self.if_gain, 0)

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq
        self.osmosdr_sink_0.set_center_freq(self.freq, 0)

    def get_bb_gain(self):
        return self.bb_gain

    def set_bb_gain(self, bb_gain):
        self.bb_gain = bb_gain
        self.osmosdr_sink_0.set_bb_gain(self.bb_gain, 0)


tb = cw_transmit()

'''
fdict holds the frequency values.
The key is the frequency at which the signal will be broadcast
Sig is a sub dictionary key which determines which signal will go out
        sig 0 is morse code (KG7WUM Call Sign)
        sig 1 is the first 11seconds of the Demo .wav file
        sig 2 is time 11-18.5 seconds
        sig 3 is time 18.5-26.0 seconds
        sig 4 is time 26.0-33.5 seconds
        sig 5 is time 33.5 to the end of the file
The frequencies are
        3.55MHz
        7.1MHz
        21.1MHz
        28.6MHz
        52MHz
        145MHz
        225MHz
        435MHz
        905MHz
slptime is a sub dictionary key which determines the time to spend at that frequency.
        These can be adjusted to fit the needs during the demo.
        
Currently the morse code sections are commented out to make it
easier to time the FM signals.
To uncomment, just delete the # ~ from the front of the line or highlight the
block of code and press ctrol + e.
'''

fdict = {
        3.55e6:{'sig':0,'slptime':22/4},      # Morse Code Signal
        7.1e6:{'sig':0,'slptime':22/4},       # Morse Code Signal
        21.1e6:{'sig':0,'slptime':22/4},      # Morse Code Signal
        28.6e6:{'sig':0,'slptime':22/4},      # Morse Code Signal
        52e6:{'sig':1,'slptime':22},
        145e6:{'sig':2,'slptime':15.0},
        225e6:{'sig':3,'slptime':15.0},
        435e6:{'sig':4,'slptime':15.0},
        905e6:{'sig':5,'slptime':15.0}
        } 
        
        
tb.start()
while True:
    for f in sorted(fdict.keys()):
        freq, sig, slptime = f, fdict[f]['sig'], fdict[f]['slptime']
        print('Setting Frequency to:', freq)
        print('Setting signal type to: ', sig)
        tb.set_signal_source(signal_source = sig)
        tb.set_freq(freq = f)
        print('Frequency has been set to: ', tb.get_freq())
        print('Signal type has been set to: ', tb.get_signal_source())
        time.sleep(slptime)
