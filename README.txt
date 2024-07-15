VIView (Last Updated: 07/14/24)
Prepared by Sankarsh Rao (srr2949@mit.edu)
NOTE: Tested on MATLAB 2021a and 2024a
If VIView is used in work, please cite our associated paper -- details included in the license.



Start-Up
=======================
To open the tool double click VIView.mlapp from the cloned Windows/MacOS/Linux directory. If this does not work, have MATLAB open and right click VIView.mlapp in the correct MATLAB directory and click "Run".

Please note that we have disabled resizing the window, since it was causing some trouble in the GUI layout.

Please wait until the logos appear before using the tool.



Default Values
=======================
The Pulse Voltage Data field has already been populated with a sample data file. This data file is one experimental voltage pulse from an APG experiment. Other sample data files are provided in the Data folder.

Checking the "Sample Gaussian" overrides the data input file and uses the Gaussian function used in Section 6 of the paper as the voltage input.



Some Things to Note
=======================

* Always start with a cable element and loads should always be connected with cable elements

* VIView forces the end of the system to be ground, i.e. V(x = L, t) = 0

* The Sample Gaussian input can be used for rapidly seeing how a waveform's reflections will look -- experimental input files may take a bit longer

* The Load Plots checkbox will plot the difference in voltage across the load and the corresponding current

* Please be patient with spark gap inputs, especially if there are multiple, as they are computationally costly to model. Please simplify the inputs if the computation takes too long (may be too many points, the max time is too long, etc.)

* The spark gap's insulating and conductive resistances are hard-coded to model air -- please feel free to change this in the back-end for other gases

* Loads have 0 length -- they are used as BCs between each cable element. As such, cable-cable interfaces have an implied load of [R,C] = [0,0]

* Checking the "Extract Solution Vecs" box before computation will save I, V, t, and x in VIViewOut.mat in the VIView directory



Example 1: Matched Load
=======================
Let's start with a simple example: a cable leading into a matched load (a resistor that matches the cable impedance exactly).

Please input the following by using the Element Type drop down menu and entering the appropriate numbers:

Cable = [10, 50, 0.6, 100]
Load = [50, 0]

The first cable is a typical 10 meter cable from your pulser to the resistor. The load is a 50 ohm resistor with no capacitance. 

Now, leave the Probe Location Field blank and make sure that the Load Plots check box is checked (to plot the properties at the load), and change the normalization parameters to be:

Char. Peak Voltage: 5e3 V
Char. Impedance: 50 Ohms
Char. Length = 10 m

Leave everything else as is, and please click Compute.

The RHS of the GUI should populate within 10-15 seconds. The Voltage vs. Time plot should be a peak at ~100 ns and and the current profile should look very similar to the voltage. Since it is a matched load, the plots should show that there are little-to-no reflections. Clicking Play under the bottom-left graph should play a video that shows this fact. Also please notice that a lot of energy is deposited in the load because it is a matched system with little-to-no reflections.

Now, let's add a probe at x = 0 (pulser) and x = 5 m (midpoint of cable). Please input 0, 5 (with the space) into the Probe Location field, do not change anything else, and click Compute again. You will see that top-left, top-right, and bottom-right plots show this change appropriately.



Example 2: Multiple Loads
=========================
Click Clear All to reset the model and GUI (this may take a couple of seconds).

Now, let's model a system with multiple loads and changing cable characteristics. Using the same voltage input file, please input the following:

Cable = [5, 75, 0.6, 100]
Cable = [2.5, 50, 0.5, 100]
Load = [500, 2e-11]
Cable = [3, 50, 0.5, 100]
Load = [250, 5e-11]
Cable = [3, 100, 0.6, 100]

Please note that the cables have different properties for no reason other than to showcase the tool's capabilities. The loads have different properties for the same reason.

Please click Compute, and once the results pop up (~20 secs since it is a complicated system) please feel free to explore the results as you wish (play the video, look at the waveforms and notice how most of the energy is deposited in the first load, how the reflection polarities differ from cable to cable, etc.).

Next, let's use the same system but instead use the sample Gaussian as an input: exp(-(t-35./33.4).^2./0.1291). This is the same sample Gaussian as used in the paper. Please click the Sample Gaussian check-box at the top, do not change anything else, and then click Compute.

Please note how this took drastically less time to produce results than the experimental case. Please also note how everything is smoother than for the experimental input (so detail is lost), but the general shape and trends are the same. As such, one can note that the Sample Gaussian is well-suited for rapidly seeing results with some loss in resolution.



Example 3: Spark Gap
=====================
Now, we can model a load that sparks/breaks down.

Click Clear All and make sure that the Sample Gaussian box is checked.

Now, input the following:

Cable = [8, 50, 0.6, 100]
Spark Gap = [5e-11, 7500] (for C, Vbreak) 

Please put 4 in the Probe Location field to see how the waveforms look in the middle of our cable. Let's also extract these solution vectors, so please click the "Extract Solution Vecs" box too.

Spark gaps are more computationally costly to model -- this run might take longer than the past cases (~30 secs) so please be patient!

Explore the waveforms and if you'd like, please also try using the saved data to get other parameters, like the cumulative current, power, etc.

This is a basic framework for understanding the tool -- please feel free to use it for your own situations. Another great use is to follow along with Case Studies 1-4, and 6 in the paper -- the outputs should match the paper exactly if the tool is used correctly!