<html>
<style>
        body { font-family: helvetica; color: black; font-size: 14pt;}
	p { padding: 0pt 5pt 0pt 10pt }
</style>
<head><title>MC-Mesh</title></head>
<body BACKGROUND="../../backgrounds/graypaper.gif">
<a href="../../index.html"><img src=../../icons/home.gif alt="Home"></a>
<a href="../index.html"><img src=../../icons/left.gif alt="Back"></a>
<hr>
<body> 
<center>
      <h1>Monte-Carlo-Based Mesh Generation</h1>
Run Java-demos for
<a href="demo/run/line.html">Line</a>, <a href="demo/run/square.html">Square</a>
<a href="demo/run/cube.html">Cube</a>, and <a href="demo/sphere/">Sphere</a>.
</center>


<h3>Summary</h3>

<p>
Mesh generation, which discretizes computation geometries into small
elements, is an important preprocessing part for numerical methods 
such as finite volume methods and finite volume methods. 
Computational domains are represented in
meshless numerical methods by nodes, or in mesh-based methods by elements
connected by nodes. Good nodal distribution is essential for the convergence of numerical methods.
In our study a new approach based on Monte Carlo simulation is proposed for
node distribution. </p>

<p>Monte Carlo (MC) methods are stochastic techniques, which use random
numbers and probability statistics to solve problems. MC methods are widely
used in various fields.
MC methods are capable of investigate complex systems such as large energy
systems, which consist of hundreds or thousands of atoms.</p>

<p>
In the proposed mesh generation method surface or volume domains to be
discretized are treated as atomic systems, and mesh nodes are considered as
interacting particles. Particles are inserted/removed into/from the system, as
well as displaced, using Grand Canonical Monte Carlo method.  With energy
minimization by the Monte Carlo simulation, a good number of particles, i.e.
nodes, are distributed with desired separation from each other. The nodal
separation is controlled by a predefined node spacing function. For mesh-based
numerical methods, well-shaped triangular mesh elements can then be generated
by connecting the nodes with constrained Delaunay methods.  </p>

<hr/>
<h3>Advantages</h3>

<p>
The GCMC method is conceptually simple and has a straightforward implementation. The resulting algorithm
 algorithm has the following advantages:
</p>

<ul>
<li>Simple robust and numerically stable.</li>
<li>Works in an identical way for both 2D and 3D cases.</li>
<li>Applicable for both mesh-based and meshless methods.</li>
</ul>

<hr/>
<h3>Comparative analysis of different algorithms</h3>

<p>
Several algorithms are implemented:
<a href="#mc">Monte Carlo (MC)</a>, 
<a href="#gcmc">Grand Canonical MC</a>,
<a href="#md">Molecular Dynamics (MD)</a>,
<a href="#md">Molecular Dynamics - Monte Carlo (MDMC)</a>,
and 
<a href="#md">Molecular Dynamics - Grand Canonical Monte Carlo (MDGCMC)</a>.
</p>

<p>The distict characteristics of the algorithms are:</p>

<ul>
<a name="mc"/>
<li>
<b>Monte-Carlo (MC):</b>
The number of nodes is fixed.
Nodes are moved with random selection of both the direction and the length of the displacement along the direction.
</li>
<a name="gcmc"/>
<li><b>Grand Cannonical Monte-Carlo (GCMC):</b>
The same as MC, but the number of nodes is selected so as to minimize the total system energy.
</li>

<a name="md"/>
<li><b>Molecular Dynamics (MD):</b>
Nodes are moved into the direction of cumulative interaction force. 
The displacements are proportional to the magnitude of the force.
</li>

<a name="mdmc"/>
<li><b>Molecular Dynamics -  Monte-Carlo (MDMC):</b> The number of nodes is fixed.
Nodes are moved into the direction of cumulative interaction force. 
The displacements along the direction are selected randomly.
</li>

<a name="mdgcmc"/>
<li><b>Molecular Dynamics - Grand Canonical Monte-Carlo (MDGCMC):</b>
The same as MDMC before, but the number of nodes is selected so as to minimize the total system energy.
</li>
</ul>

<p>
As can be seen from demo examples below, MD method is quite impractical for node placement. MC method leads to a more accurate and uniform node placement, compared to MDMC. However the latter is much faster. Likewise, MDGCMC method is fater than GCMC, but can lead to void regions. The optimum strategy can be running the MDGCMC algorithm at the beginning and GCMC in the final stages.
</p>


<hr/>
<h3>Demonstration examples</h3>

<p>
In the demo examples the Lennart Jones
potential is used to control the dynamics of the particles.
</p>
<p>
Java-applets with the demonstration of the method can be run for
<a href="demo/run/line.html">Line</a>, <a href="demo/run/square.html">Square</a>
<a href="demo/run/cube.html">Cube</a>, and <a href="demo/sphere/">Sphere</a>.
</p>

<h4>Source codes</h4>

<p>Here are the <a href="demo/src">the source codes of the demo examples,</a> which are also available as a 
<a href="../../download/mcmesh.tgz">tar-gzipped file</a>.
</p>

<hr/>
<h3>Bibliography</h3>

<?php

$user = 'andrei';
$password = 'zaq12wsx';
$database = 'mulphys';

$link = mysql_connect('localhost', $user, $password);
@mysql_select_db($database) or die( "Unable to select database");

$query="select * from bib where author like '%Smirnov%' && (title like '%mesh%' || title like '%node%') order by year desc";

	$result = mysql_query($query);
	if (!$result)
	{	die ("Query failed");
	}
	$n = mysql_numrows($result);
	for ($i = 0; $i < $n; $i++)
	{
		$title = mysql_result($result, $i, "title");
		if($title!='') echo "<b>$title</b><BR>";
		$author = mysql_result($result, $i, "author");
		if($author!='') echo "$author<br>";
		$type = mysql_result($result, $i, "type");
		if(($selected_type=='any'||$selected_type=='other')&&$type!='') echo "Document type: $type<br>";
		$year = mysql_result($result, $i, "year");
		if($year!='') echo "Year: <font color='brown'><b>$year</b></font><BR>";
		$booktitle = mysql_result($result, $i, "booktitle");
		if($booktitle!='') echo "Conference: $booktitle<BR>";
		$journal = mysql_result($result, $i, "journal");
		if($journal!='') echo "Journal: $journal<BR>";
		$address = mysql_result($result, $i, "address");
		if($address!='') echo "Address: $address<BR>";
		$series = mysql_result($result, $i, "series");
		if($series!='') echo "Series: $series<BR>";
		$editor = mysql_result($result, $i, "editor");
		if($editor!='') echo "Editor: $editor<BR>";
		$isbn = mysql_result($result, $i, "isbn");
		if($isbn!='') echo "ISBN: $isbn<br>";
		$volume = mysql_result($result, $i, "volume");
		if($volume!='') echo "Volume: $volume<br>";
		$number = mysql_result($result, $i, "number");
		if($number!='') echo "Number: $number<BR>";
		$publisher = mysql_result($result, $i, "publisher");
		if($publisher!='') echo "Publisher: $publisher<BR>";
		$pages = mysql_result($result, $i, "pages");
		if($pages!='') echo "Pages: $pages<BR>";
		$url = mysql_result($result, $i, "url");
		if($url!='') echo "URL: <a href=\"$url\">$url</a><br>";
		echo "<hr>";
	}

?>

</body>
</html>	  

