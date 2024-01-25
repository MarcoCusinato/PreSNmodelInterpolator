# PreSNmodelInterpolator
## Introduction
Small `Python` package to interpolate pre-supernova files from [MESA](https://docs.mesastar.org/en/release-r23.05.1/#) or [KEPLER](https://2sn.org/kepler/doc/index.html) to Aenus format. The package will generate seven files:
 - `Heger.pars`: Fortran namelist with the number of radius cells;
 - `initial_model.x.dat`:  4 column, cell number, left cell radius, central cell radius and right cell radius;
 - `initial_model.y.dat`: two columns and one lines, representing the angular dimension;
 - `nuclei.dat`: matter composition of the model
 - `nuclei.pars`: Fortran namelist containing the information on the matter composition;
 - `star.dat`: thermodynamics variables of the model;
 - `star.pars`: model information, such as maximum and minimum radius, number of grid cells and enclosed mass.
## Requirements
The following `Python` package is mandatory to run the script:
 - [NumPy](https://numpy.org/)
## How to run
Before starting to use the package make sure that the original presupernova model is made of just one file. Running
```
python generate_presn_model.py --model-path /path/to/model
```
will produce the seven files mentioned above with 16000 radial cells, from 0 to $10^{13}$ cm in logscale, with the first cell being $8\cdot10^{4}$. To modify the parameters please run
```
python generate_presn_model.py -h
```
## What's in the files?
### star.dat
The thermodynamics quantities appearing in the `star.dat` are, in order: <br>
<table>
  <tr>
    <td colspan="11">Quantities present in every model </td>
    <td colspan="2">Pesent only if the model has a magnetic field</td>
  </tr>
  <tr>
    <td>Cell number</td>
    <td>Radius</td>
    <td>Density</td>
    <td>Temperature</td>
    <td>$Y_e$</td>
    <td>Pressure</td>
    <td>Entropy</td>
    <td>Internal energy</td>
    <td>$\overline{A}$</td>
    <td>Radial velocity</td>
    <td>Angular velocity</td>
    <td>Toroidal B field</td>
    <td>Poloidal B field</td>
  </tr>
</table>

## nuclei.dat

The composition quantities appearing in the composition file are the 20 representative species of nuclei. 
<table>
  <tbody><tr>
    <th rowspan="2">column</th>
    <th colspan="2">network</th>
  </tr>
  <tr>
    <th> APPROX </th>
    <th> QSE / NSE</th>
  </tr>
  <tr>
    <td><code>n</code></td>
    <td>cell number</td>
    <td>cell number</td>
  </tr>
  <tr>
    <td><code>r</code></td>
    <td>radius</td>
    <td>radius</td>
  </tr>
  <tr>
    <td><code>neutrons</code></td>
    <td>neutrons</td>
    <td>neutrons</td>
    
  </tr>
  <tr>
    <td><code>H1</code></td>
    <td><sup>1</sup>H</td>
    <td>protons</td>
    
  </tr>
  <tr>
    <td><code>He3</code></td>
    <td><sup>3</sup>He</td>
    <td>---</td>
    
  </tr>
  <tr>
    <td><code>He4</code></td>
    <td><sup>4</sup>He</td>
    <td>1 &lt; A &lt; 6</td>
    
  </tr>
  <tr>
    <td><code>C12</code></td>
    <td><sup>12</sup>C</td>
    <td>---</td>
    
  </tr>
  <tr>
    <td><code>N14</code></td>
    <td><sup>14</sup>N</td>
    <td>---</td>
    
  </tr>
  <tr>
    <td><code>O16</code></td>
    <td><sup>16</sup>O</td>
    <td><sup>16</sup>O (unburned yet)</td>
    
  </tr>
  <tr>
    <td><code>Ne20</code></td>
    <td><sup>20</sup>Ne</td>
    <td>---</td>
    
  </tr>
  <tr>
    <td><code>Mg24</code></td>
    <td><sup>24</sup>Mg</td>
    <td>22 &lt; A &lt; 29, excluding <sup>28</sup>Si</td>
    
  </tr>
  <tr>
    <td><code>Si28</code></td>
    <td><sup>28</sup>Si</td>
    <td><sup>28</sup>Si</td>
    
  </tr>
  <tr>
    <td><code>S32</code></td>
    <td><sup>32</sup>S</td>
    <td>28 &lt; A &lt; 36</td>
    
  </tr>
  <tr>
    <td><code>Ar36</code></td>
    <td><sup>36</sup>Ar</td>
    <td>35 &lt; A &lt; 40</td>
    
  </tr>
  <tr>
    <td><code>Ca40</code></td>
    <td><sup>40</sup>Ca</td>
    <td>39 &lt; A &lt; 44</td>
    
  </tr>
  <tr>
    <td><code>Ti44</code></td>
    <td><sup>44</sup>Ti</td>
    <td>43 &lt; A &lt; 48</td>
    
  </tr>
  <tr>
    <td><code>Cr48</code></td>
    <td><sup>48</sup>Cr</td>
    <td>47 &lt; A &lt; 52</td>
    
  </tr>
  <tr>
    <td><code>Fe52</code></td>
    <td><sup>52</sup>Fe</td>
    <td>---</td>
    
  </tr>
  <tr>
    <td><code>Fe54</code></td>
    <td><sup>54</sup>Fe (+<sup>56</sup>Fe)</td>
    <td>A approximately 2*Z+2, Iron Peak</td>
    
  </tr>
  <tr>
    <td><code>Ni56</code></td>
    <td><sup>56</sup>Ni</td>
    <td>A &lt; 2*Z+2, Iron Peak</td>
    
  </tr>
  <tr>
    <td><code>Fe56</code></td>
    <td>---</td>
    <td><sup>56</sup>Fe only</td>
    
  </tr>
  <tr>
    <td><code>'Fe'</code></td>
    <td>---</td>
    <td>A &gt; 2*Z + 3, Iron Peak, excluding <sup>56</sup>Fe</td>
    
  </tr>
</tbody></table>
