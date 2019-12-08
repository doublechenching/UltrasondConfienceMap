General
=======

C++ code for confidence estimation using random walks.
For questions and feedback please contact Athanasios Karamalis (karamali@in.tum.de).

Ultrasound confidence maps
A. Karamalis, W. Wein, T. Klein, N. Navab: Ultrasound Confidence Maps
using Random Walks, Medical Image Analysis, 16, 6, 1101 - 1112, 2012
DOI: http://dx.doi.org/10.1016/j.media.2012.07.005

Chair for Computer Aided Medical Procedures (CAMP)
Technische Universität München
Written by: Athanasios Karamalis
Email: karamali@in.tum.de

THE WORK IS FOR RESEARCH AND NON-COMMERCIAL PURPOSES ONLY.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL CAMP BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


The code further implements the random walks for image segmentation algorithm from:
Grady, L.: Random walks for image segmentation.IEEE Transactions on Pattern Analysis and Machine Intelligence 28(2006), 1768–1783.

Building code
=============

Use CMake to create platform specific build with the compiler and IDE of your choice.

Code requires Eigen eigen3.1rc1 or above and ViennaCL-1.3.0 or above.

Alternatively, you can of course remove the ViennaCL GPU solver. This solver is intended for the 3D random walks segmentation problem (also implemented in this code).
For 2D problems use direct solvers like LLT (see below).

Usage
=====

Use the ConfidenceMaps2DFacade, which is a facade to the underlying confidence estimation code. 
To solve the system of linear equations use the direct Eigen-LLT. Alternatively, in the unlikely case that you
run out of memory use the iterative solver from SparseSolverEigenBiCGSTAB.

Attention must be paid when providing the input images. The code assumes images to be provided in a linear array (std::vector) with COLUMN-MAJOR order.
The same format is used to provide the confidence map output.


LICENSE
=======

Code uses the Eigen library. Check LICENSE file for more information.

