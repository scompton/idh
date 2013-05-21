package util;

// Translated from C++ Version:
// nvModel.h - Model support class
//
// The nvModel class implements an interface for a multipurpose model
// object. This class is useful for loading and formatting meshes
// for use by OpenGL. It can compute face normals, tangents, and
// adjacency information. The class supports the obj file format.
//
// Author: Evan Hart
// Email: sdkfeedback@nvidia.com
//
// Copyright (c) NVIDIA Corporation. All rights reserved.
////////////////////////////////////////////////////////////////////////////////

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Vector;

public class Model {

	public Model CreateModel() {
		return new Model();
	}

	public static Model load(String modelFilename) {
    Model model = new Model();
    System.err.println("loading OBJ...\n");
    model.loadModelFromFile( modelFilename );
    System.err.println("compiling mesh...\n");
    model.compileModel();
    System.err.println(model.getPositionCount() + " vertices");
    System.err.println((model.getIndexCount()/3) + " triangles");
    return model;
  }
	
	public class IdxSet {
		Integer pIndex = 0;
		Integer nIndex = 0;

		boolean lessThan (IdxSet rhs) {
		  if (pIndex < rhs.pIndex)
				return true;
			else if (pIndex == rhs.pIndex) {
				if (nIndex < rhs.nIndex)
					return true;
			}
			return false;
		}
	};

	//
	// loadModelFromFile
	//
	//    This function attempts to determine the type of
	//  the filename passed as a parameter. If it understands
	//  that file type, it attempts to parse and load the file
	//  into its raw data structures. If the file type is
	//  recognized and successfully parsed, the function returns
	//  true, otherwise it returns false.
	//
	//////////////////////////////////////////////////////////////
	public boolean loadModelFromFile( String file ) {
	  BufferedReader input = null;
	  try {
	    input = new BufferedReader(
	        new InputStreamReader(new FileInputStream(file)));
	    String line = null;
	    float[] val = new float[4];
	    int[][] idx = new int[3][3];
	    boolean hasNormals = false;

	    while ( (line = input.readLine()) != null) {
	      if (line.isEmpty())
	        continue;
	      switch (line.charAt(0)) {
	      case '#':
	        break;

	      case 'v':
	        switch (line.charAt(1)) {

	        case ' ':
	          line = line.substring(line.indexOf(" ") + 1).trim();
	          //vertex, 3 or 4 components
	          val[0] = Float.valueOf(line.substring(0,line.indexOf(" ")));
	          line = line.substring(line.indexOf(" ") + 1).trim();
	          val[1] = Float.valueOf(line.substring(0,line.indexOf(" " )));
	          line = line.substring(line.indexOf(" ") + 1 ).trim();
	          val[2] = Float.valueOf(line);
	          _positions.add(val[0]);
	          _positions.add(val[1]);
	          _positions.add(val[2]);
	          break;

	        case 'n':
	          //normal, 3 components
	          line = line.substring(line.indexOf(" ") + 1).trim();
	          val[0] = Float.valueOf(line.substring(0,line.indexOf( " " )));
	          line = line.substring(line.indexOf(" ") + 1).trim();
	          val[1] = Float.valueOf(line.substring(0,line.indexOf( " " )));
	          line = line.substring(line.indexOf(" ") + 1).trim();
	          val[2] = Float.valueOf(line);
	          _normals.add(val[0]);
	          _normals.add(val[1]);
	          _normals.add(val[2]);
	          break;
	        
	        default:
            break;
	        }
	        break;
	        

	      case 'f':
	        //face
	        line = line.substring(line.indexOf(" ") + 2);

	        idx[0][0] = Integer.valueOf(
	            line.substring(0,line.indexOf("//"))).intValue();
	        line = line.substring(line.indexOf("//") + 2);
	        idx[0][1] = Integer.valueOf(
	            line.substring(0,line.indexOf(" "))).intValue();

	        //This face has vertex and normal indices
	        // in .obj, f v1 .... the vertex index used start from 1, so -1 here
	        //remap them to the right spot
	        idx[0][0] = (idx[0][0] > 0) ?
	            (idx[0][0] - 1) : (_positions.size() - idx[0][0]);
	        idx[0][1] = (idx[0][1] > 0) ?
	            (idx[0][1] - 1) : (_normals.size() - idx[0][1]);

	        //grab the second vertex to prime
	        line = line.substring(line.indexOf(" ") + 1);
	        idx[1][0] = Integer.valueOf(line.substring(0,line.indexOf("//")));
	        line = line.substring(line.indexOf("//") + 2);
	        idx[1][1] = Integer.valueOf(line.substring(0,line.indexOf(" ")));

	        //remap them to the right spot
	        idx[1][0] = (idx[1][0] > 0) ? 
	            (idx[1][0] - 1) : (_positions.size() - idx[1][0]);
	        idx[1][1] = (idx[1][1] > 0) ? 
	            (idx[1][1] - 1) : (_normals.size() - idx[1][1]);

	        //grab the third vertex to prime
	        line = line.substring(line.indexOf(" ") + 1);
	        idx[2][0] = Integer.valueOf(line.substring(0,line.indexOf("//")));
	        line = line.substring(line.indexOf("//") + 2);
	        idx[2][1] = Integer.valueOf(line);
	          
	        //remap them to the right spot
	        idx[2][0] = (idx[2][0] > 0) ? 
	            (idx[2][0] - 1) : (_positions.size() - idx[2][0]);
	        idx[2][1] = (idx[2][1] > 0) ? 
	            (idx[2][1] - 1) : (_normals.size() - idx[2][1]);

	        //add the indices
	        for (int ii = 0; ii < 3; ii++) {
	          _pIndex.add( idx[ii][0]);
	          _nIndex.add( idx[ii][1]);
	        }
	        //prepare for the next iteration, the num 0 does not change.
	        idx[1][0] = idx[2][0];
	        idx[1][1] = idx[2][1];
	        hasNormals = true;
	        break;

	      default:
	        break;
	      };								
	    }			
	    //post-process data
	    //free anything that ended up being unused
	    if (!hasNormals) {
	      _normals.clear();
	      _nIndex.clear();
	    }

	    _posSize = 3;
	    return true;

	  } catch (FileNotFoundException kFNF) {
	    System.err.println("Unable to find the shader file " + file);
	  } catch (IOException kIO) {
	    System.err.println("Problem reading the shader file " + file);
	  } catch (NumberFormatException kIO) {
	    System.err.println("Problem reading the shader file " + file);
	  } finally {
	    try {
	      if (input != null) {
	        input.close();
	      }
	    } catch (IOException closee) {}
	  }
		return false;
	}

	//
	//  compileModel
	//
	//    This function takes the raw model data in the internal
	//  structures, and attempts to bring it to a format directly
	//  accepted for vertex array style rendering. This means that
	//  a unique compiled vertex will exist for each unique
	//  combination of position, normal, tex coords, etc that are
	//  used in the model. The prim parameter, tells the model
	//  what type of index list to compile. By default it compiles
	//  a simple triangle mesh with no connectivity. 
	//
	public void compileModel( ) {
		_v = new float[_pIndex.size()*3];
		_n = new float[_pIndex.size()*3];
		for (int i=0; i<_pIndex.size(); i++) {
			IdxSet idx = new IdxSet();
			idx.pIndex = _pIndex.elementAt(i);
			if ( _normals.size() > 0)
				idx.nIndex = _nIndex.elementAt(i);
			else
				idx.nIndex = 0;

			for (int j=0; j<3; j++) {
			  _v[i*3+j] = _positions.elementAt(idx.pIndex*_posSize+j);
			  _n[i*3+j] = _normals.elementAt(idx.nIndex*3+j);
			}
		}
	}

	public int getPositionCount()  {
		return (_posSize > 0) ? _positions.size() / _posSize : 0;
	}

	public int getNormalCount()  {
		return _normals.size()/3;
	}

	public int getIndexCount()  {
		return _pIndex.size();
	}

	public float[] getVertices() {
    return _v;
  }

	public float[] getNormals() {
    return _n;
  }
	
	private Vector<Float> _positions = new Vector<Float>();
	private Vector<Float> _normals = new Vector<Float>();
	private Vector<Integer> _pIndex = new Vector<Integer>();
	private Vector<Integer> _nIndex = new Vector<Integer>();

	private float[] _v;
	private float[] _n;
	private int _posSize;
};
