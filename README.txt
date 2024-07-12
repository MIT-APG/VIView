VIView (Last Updated: 06/20/24)
Prepared by Sankarsh Rao (srr2949@mit.edu)
NOTE: For MATLAB2024a onward! Might run into errors with previous versions



Start-Up
=======================
To open the tool, if you do not already have MATLAB open, double click VIView.mlapp. If this does not work, have MATLAB open and right click VIView.mlapp in the directory window and click "Run".

Please wait until the logos appear before using the tool.



Default Values
=======================
The Pulse Voltage Data field has already been populated with a sample data file. This data file is an experimental voltage waveform from an APG experiment provided in Figure 3 in the paper. Other sample data files are provided in the Data folder.

The scaling parameters at the bottom of the GUI are also filled in with reasonable values for the simulation.

Checking the "Sample Gaussian" overrides the data input file and uses the Gaussian function used in Section 4 of the paper as the voltage input.



Best-Practices
=======================
Some things to note:

* Always start with a cable element

* Always end with a cable element -- this is how the ground end condition is applied. If your system has a grounding cable already, great! If not, we suggest you put in a small 1cm cable at the end of the system, as it will not noticeably change the resulting waveforms: Cable = [1e-2, 50, 0.6, 1]

* The Sample Gaussian input can be used for rapidly seeing how a waveform's reflections will look -- experimental input files may take longer

* The Load Plots checkbox will plot the difference in voltage across the load and the corresponding current

* Depending on the input the model may take ~minutes, especially for spark gap inputs -- please be patient! If it takes too long (i.e. hangs at 2 or 3/5 in the progress bar, or gives an error message above the progress bar), please do control + C to end the processing, re-open VIView and please simplify the inputs. There may be too many points, the max time is too long, etc.

* All loads are constant RC loads and the spark gap voltage threshold is hard-coded to model air -- please feel free to change this in the back-end for other gases

* Cable inputs should have non-zero inputs otherwise singularities will occur

* Loads have 0 length -- they are used as BCs between each cable element. As such, cable-cable interfaces have an implied load of [R,C] = [0,0]

* Checking the "Extract Solution Vecs" box before computation will save I, V, t, and x as a VIViewOut.mat in the VIView directory



Example 1: Matched Load
=======================
Let's start with a simple example: a cable leading into a matched load (a resistor that matches the cable impedance exactly), leading to ground.

Please keep the default values and please input the following by using the Element Type drop down menu and entering the appropriate numbers:

Cable = [5, 50, 0.6, 200]
Load = [50, 0]
Cable = [1e-2, 50, 0.6, 1]

The first cable is a typical 5 meter cable from your power supply to the resistor. The load is a 50 ohm resistor with no capacitance. As all systems for the model need to start and end with a cable, a very small cable with the same properties was included at the end, leading to ground.

Now, check the box for Load Plots (to plot the properties at the load), do not change any of the scaling parameters and click Compute.

The RHS of the GUI should populate within a minute. The Voltage vs. Time plot should be a peak at ~7e-8 secs and and the current profile should look very similar to the voltage. Since it is a matched load, the plots should show that there are little-to-no reflections. Clicking Play under the bottom-left graph should play a video that shows this fact. Please wait until the video is done playing before inputting any other commands.

Now, let's add a probe at x = 0 and x = 2.5 m. Please input 0, 2.5 (with the space) into the Probe Location field, do not change anything else, and click Compute again. You will see that top-left, top-right, and bottom-right plots show this change appropriately.



Example 2: Multiple Loads
=========================
Click Clear All to reset the model and GUI.

Now, let's model a system with multiple loads and changing cable characteristics. Using the same voltage input file, please input the following:

Cable = [5, 75, 0.6, 100]
Cable = [2.5, 50, 0.5, 100]
Load = [5000, 2e-11]
Cable = [3, 50, 0.5, 100]
Load = [2500, 5e-11]
Cable = [3, 100, 0.6, 100]

Sparking is not modeled in this case. Also, please notice how the grounding cable is long for no reason other than to showcase the tool's capabilities. The loads have different properties for the same reason.

Check the Load Plots button once again, click Compute, and once the results pop up please feel free to explore the results as you wish (play the video, look at the waveforms and notice how most of the energy is deposited in the first load).

Next, let's use the same system but instead use the sample Gaussian as an input: (exp(-(t-1.5).^2/0.2)). This is the same sample Gaussian as used in the paper. Please click the Sample Gaussian check-box, and do not change anything else, and then click Compute.

Please note how this took less time to produce results than the experimental case. Please also look at the plots and note how everything is smoother than for the experimental input (so detail is lost), but the general shape and trends are the same. As such, one can note that the Sample Gaussian is well-suited for rapidly seeing results with some loss in resolution.



Example 3: Spark Gap
=====================
Now, we can model a load that sparks/breaks down.

Click Clear All and uncheck the Sample Gaussian box so that we use the default experimental file.

Now, input the following:

Cable = [4, 50, 0.6, 100]
Spark Gap = [default]
Cable = [1e-2, 50, 0.6, 1] 

Let's see what the waveforms look like at say, 3 m so please put 3 in the Probe Location field. Also, please check the Load Plots checkbox, and change the Max Time parameter to 400e-09 to shorten the max time so that the simulation does not take too long. Spark gaps are more computationally costly to model, so this run might take a couple of minutes to complete. It might hang on 2/5 in the progress bar, but please be patient!

Let's also extract these solution vectors, so please click the "Extract Solution Vecs" box too.

Now, click Compute. BUT WAIT, WE HIT AN ERROR! The message above the progress bar describes the error -- when you have an ODE solve error, please vary the # of points in the cables. This is the price we pay for generalizing systems and using MATLAB's ode15s, which sometimes does not like the resolution we give it. It sometimes throws errors if you put in 1 point vs. 2 points, so we encourage the user to vary the # of points if an error is reached.

To fix the error, let's change the # of points in the grounding cable to 2 by first clicking Delete Element and adding this element:

Cable = [1e-2, 50, 0.6, 2] 

Keep everything else the same and click Compute and you should see the results on the RHS of the GUI!

Explore the waveforms and if you'd like, please also try using the saved data to get other parameters, like the cumulative current, power, etc.

This is a basic framework for understanding the tool -- please feel free to use it for your own situations. Another great use is to follow along with Case Studies 1-4, and 6 in the paper -- the outputs should match the paper exactly if the tool is used correctly!