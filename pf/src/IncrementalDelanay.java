/********* Incremental Delaunay triangulation begins here            *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  boundingbox()   Form an "infinite" bounding triangle to insert vertices  */
/*                  into.                                                    */
/*                                                                           */
/*  The vertices at "infinity" are assigned finite coordinates, which are    */
/*  used by the point location routines, but (mostly) ignored by the         */
/*  Delaunay edge flip routines.                                             */
/*                                                                           */
/*****************************************************************************/

void boundingbox(mesh m) //, struct behavior *b)
{
//-	otri inftri=new otri();          /* Handle for the triangular bounding box. */
	REAL width;

//-  if (b->verbose) {
System.out.println("  Creating triangular bounding box.");
//-  }
  /* Find the width (or height, whichever is larger) of the triangulation. */
  width = m.xmax - m.xmin;
  if (m.ymax - m.ymin > width) {
    width = m.ymax - m.ymin;
  }
  if (width == 0.0F) {
    width = 1.0F;
  }
  /* Create the vertices of the bounding box. */
//-  m.infvertex1 = (vertex) trimalloc(m.vertices.itembytes);
//-  m.infvertex2 = (vertex) trimalloc(m.vertices.itembytes);
//-  m.infvertex3 = (vertex) trimalloc(m.vertices.itembytes);
  m.infvertex1.x[0] = m.xmin - 50.0F * width;
  m.infvertex1.x[1] = m.ymin - 40.0F * width;
  m.infvertex2.x[0] = m.xmax + 50.0F * width;
  m.infvertex2.x[1] = m.ymin - 40.0F * width;
  m.infvertex3.x[0] = 0.5F * (m.xmin + m.xmax);
  m.infvertex3.x[1] = m.ymax + 60.0F * width;

  /* Create the bounding box. */
  //-maketriangle(m, b, &inftri);
	otri inftri=maketriangle(m);
  setorg(inftri, m.infvertex1);
  setdest(inftri, m.infvertex2);
  setapex(inftri, m.infvertex3);
  /* Link dummytri to the bounding box so we can always find an */
  /*   edge to begin searching (point location) from.           */
//-  m.dummytri[0] = inftri.tri;
  m.dummytri = inftri.tri;
//- if (b->verbose > 2) {
System.out.println("  Creating ");
//-    printtriangle(m, b, &inftri);
//-  }
}

/*****************************************************************************/
/*                                                                           */
/*  removebox()   Remove the "infinite" bounding triangle, setting boundary  */
/*                markers as appropriate.                                    */
/*                                                                           */
/*  The triangular bounding box has three boundary triangles (one for each   */
/*  side of the bounding box), and a bunch of triangles fanning out from     */
/*  the three bounding box vertices (one triangle for each edge of the       */
/*  convex hull of the inner mesh).  This routine removes these triangles.   */
/*                                                                           */
/*  Returns the number of edges on the convex hull of the triangulation.     */
/*                                                                           */
/*****************************************************************************/
long removebox(mesh m) //, behavior *b)
{
otri deadtriangle=new otri(),
	searchedge=new otri(),
	checkedge=new otri(),
	nextedge=new otri(), finaledge=new otri(), dissolveedge=new otri();
  vertex markorg=new vertex();
  long hullsize;
  triangle ptr=new triangle();                         /* Temporary variable used by sym(). */

//-  if (b->verbose) {
System.out.println("  Removing triangular bounding box.");
//-  }
  /* Find a boundary triangle. */
  nextedge.tri = m.dummytri;
  nextedge.orient = 0;
  symself(nextedge);
  /* Mark a place to stop. */
  lprev(nextedge, finaledge);
  lnextself(nextedge);
  symself(nextedge);
  /* Find a triangle (on the boundary of the vertex set) that isn't */
  /*   a bounding box triangle.                                     */
  lprev(nextedge, searchedge);
  symself(searchedge);
  /* Check whether nextedge is another boundary triangle */
  /*   adjacent to the first one.                        */
  lnext(nextedge, checkedge);
  symself(checkedge);
  if (checkedge.tri == m.dummytri) {
    /* Go on to the next triangle.  There are only three boundary   */
    /*   triangles, and this next triangle cannot be the third one, */
    /*   so it's safe to stop here.                                 */
    lprevself(searchedge);
    symself(searchedge);
  }
  /* Find a new boundary edge to search from, as the current search */
  /*   edge lies on a bounding box triangle and will be deleted.    */
//-	m.dummytri[0] = encode(searchedge);
  m.dummytri = encode(searchedge).tri;
  hullsize = -2l;
  while (!otriequal(nextedge, finaledge)) {
    hullsize++;
    lprev(nextedge, dissolveedge);
    symself(dissolveedge);
    /* If not using a PSLG, the vertices should be marked now. */
    /*   (If using a PSLG, markhull() will do the job.)        */
//-    if (!b->poly) {//poly=0 by default
      /* Be careful!  One must check for the case where all the input     */
      /*   vertices are collinear, and thus all the triangles are part of */
      /*   the bounding box.  Otherwise, the setvertexmark() call below   */
      /*   will cause a bad pointer reference.                            */
      if (dissolveedge.tri != m.dummytri) {
        org(dissolveedge, markorg);
        if (vertexmark(markorg) == 0) {
          setvertexmark(markorg, 1);
        }
      }
//-    }
    /* Disconnect the bounding box triangle from the mesh triangle. */
    dissolve(dissolveedge);
    lnext(nextedge, deadtriangle);
    sym(deadtriangle, nextedge);
    /* Get rid of the bounding box triangle. */
    triangledealloc(m, deadtriangle.tri);
    /* Do we need to turn the corner? */
    if (nextedge.tri == m.dummytri) {
      /* Turn the corner. */
      otricopy(dissolveedge, nextedge);
    }
  }
  triangledealloc(m, finaledge.tri);

//-  trifree((VOID *) m->infvertex1);  /* Deallocate the bounding box vertices. */
//-  trifree((VOID *) m->infvertex2);
//-  trifree((VOID *) m->infvertex3);
	m.infvertex1=null;
	m.infvertex2=null;
	m.infvertex3=null;

  return hullsize;
}
/*****************************************************************************/
/*                                                                           */
/*  incrementaldelaunay()   Form a Delaunay triangulation by incrementally   */
/*                          inserting vertices.                              */
/*                                                                           */
/*  Returns the number of edges on the convex hull of the triangulation.     */
/*                                                                           */
/*****************************************************************************/
long incrementaldelaunay(mesh m) //, struct behavior *b)
{
	otri starttri=new otri();
//-  vertex vertexloop; //=new vertex();

  /* Create a triangular bounding box. */
  boundingbox(m);
//-  if (b->verbose) {
System.out.println("  Incrementally inserting vertices.");
//-  }
//-  traversalinit(m.vertices);
//-  vertexloop = vertextraverse(m);
	m.vertices.goFirst();
//-  while (vertexloop != (vertex) NULL) {
	do {
		starttri.tri = m.dummytri;
		vertex vertexloop=m.vertices.getItem();
//-    if (insertvertex(m, b, vertexloop, &starttri, (struct osub *) NULL, 0, 0)
    if (insertvertex(m, vertexloop, starttri, null, 0, 0)
        == insertvertexresult.DUPLICATEVERTEX) {
//-      if (!b->quiet) {
System.out.println("Warning:  A duplicate vertex at (%.12g, %.12g) appeared and was ignored.");
//-               vertexloop[0], vertexloop[1]);
//-      }
      setvertextype(vertexloop, VertexType.UNDEADVERTEX.index());
      m.undeads++;
    }
//-    vertexloop = vertextraverse(m);
		m.vertices.goNext();
	}	while(!m.vertices.isFirst());
  /* Remove the bounding box. */
  return removebox(m);
}
/**                                                                         **/
/**                                                                         **/
/********* Incremental Delaunay triangulation ends here              *********/

/********* Sweepline Delaunay triangulation begins here              *********/
/**                                                                         **/
/**                                                                         **/
void eventheapinsert(event[] heap, int heapsize, event newevent)
{
  REAL eventx, eventy;
  int eventnum;
  int parent;
  boolean notdone;

  eventx = newevent.xkey;
  eventy = newevent.ykey;
  eventnum = heapsize;
  notdone = eventnum > 0;
  while (notdone) {
    parent = (eventnum - 1) >> 1;
    if ((heap[parent].ykey < eventy) ||
        ((heap[parent].ykey == eventy)
         && (heap[parent].xkey <= eventx))) {
      notdone = false;
    } else {
      heap[eventnum] = heap[parent];
      heap[eventnum].heapposition = eventnum;

      eventnum = parent;
      notdone = eventnum > 0;
    }
  }
  heap[eventnum] = newevent;
  newevent.heapposition = eventnum;
}
void eventheapify(event[] heap, int heapsize, int eventnum)
{
event thisevent;
REAL eventx, eventy;
int leftchild, rightchild;
int smallest;
boolean notdone;

  thisevent = heap[eventnum];
  eventx = thisevent.xkey;
  eventy = thisevent.ykey;
  leftchild = 2 * eventnum + 1;
  notdone = leftchild < heapsize;
  while (notdone) {
    if ((heap[leftchild].ykey < eventy) ||
        ((heap[leftchild].ykey == eventy)
         && (heap[leftchild].xkey < eventx))) {
      smallest = leftchild;
    } else {
      smallest = eventnum;
    }
    rightchild = leftchild + 1;
    if (rightchild < heapsize) {
      if ((heap[rightchild].ykey < heap[smallest].ykey) ||
          ((heap[rightchild].ykey == heap[smallest].ykey)
           && (heap[rightchild].xkey < heap[smallest].xkey))) {
        smallest = rightchild;
      }
    }
    if (smallest == eventnum) {
      notdone = false;
    } else {
      heap[eventnum] = heap[smallest];
      heap[eventnum].heapposition = eventnum;
      heap[smallest] = thisevent;
      thisevent.heapposition = smallest;

      eventnum = smallest;
      leftchild = 2 * eventnum + 1;
      notdone = leftchild < heapsize;
    }
  }
}
void eventheapdelete(event[] heap, int heapsize, int eventnum)
{
event moveevent;
REAL eventx, eventy;
int parent;
boolean notdone;

  moveevent = heap[heapsize - 1];
  if (eventnum > 0) {
    eventx = moveevent.xkey;
    eventy = moveevent.ykey;
    do {
      parent = (eventnum - 1) >> 1;
      if ((heap[parent].ykey < eventy) ||
          ((heap[parent].ykey == eventy)
           && (heap[parent].xkey <= eventx))) {
        notdone = false;
      } else {
        heap[eventnum] = heap[parent];
        heap[eventnum].heapposition = eventnum;

        eventnum = parent;
        notdone = eventnum > 0;
      }
    } while (notdone);
  }
  heap[eventnum] = moveevent;
  moveevent.heapposition = eventnum;
  eventheapify(heap, heapsize - 1, eventnum);
}






