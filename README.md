# Zemax-Spherical-Thermal
First order thermal modeling of spherical GRIN

NOTE (4/3/2019): Sincere apologies to early adopters.  It was discovered just today that 
the original upload (7/17/2018) was inadvertently corrupted by an intermediate version of
the code.  The thermal variation of geometric parameters was modeled, but thermal variation
of material refractive indices had yet to be implemented. The current documents and DLL
represent what was *intended* to be uploaded.  Again, apologies for the confusion.

Bare-bones release July 17, 2018.  This represents an extension of the spherical GRIN model
available at https://github.com/GrinLens/Zemax-Spherical-GRIN/releases.  It is highly
recommended that one read the documentation available there before attempting to use this tool.
Example files may follow.

It should be acknowledged that the thermal treatment of any lens in Zemax is awkward
and prone to error.  This is doubly so when handling GRIN lenses.  It is recommended that,
at the very least, an initial user of this tool model a GRIN-based singlet with a homogeneous
composition and compare it to an identical, all-Zemax homogeneous singlet.  Ensure that ray
traces and lens specifications of the two lenses agree, in detail, at multiple temperatures
before embarking on GRIN design.

NOTE: This first-order model should be treated as a guide to order-of-magnitude thermal trends,
only. The authors can make no claim to rigorous accuracy for fielded applications, which requires
detailed, material-specific information not provided in this general package.

# Acknowledgements
This work was developed at the Naval Research Laboratory by (in alphabetical order): Guy
Beadie, Richard Flynn, Casey Kretzer, and Armand Rosenberg. Many others contributed via
their valuable feedback. In particular, we would like to acknowledge the input of Michael
Brindza, Ozan Cakmakci, Howard Fein, Erin Fleet, Predrag Milojkovic, Michael Ponting, John
Rogers, and Yehudi Self-Medlin. The principal funding was provided by the Office of Naval
Research, with initial support from DARPA’s Manufacturable Gradient Index Optics (M-GRIN)
program.

# Disclaimer
2018, US Naval Research Laboratory (NRL)
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
in compliance with the License. You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License
is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the specific language governing
permissions and limitations under the License.
