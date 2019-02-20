#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    
    int length = width * height;
    state.image_color = new pixel[length];

    for(int i = 0; i < length; i++){
	state.image_color[i] = make_pixel(0, 0, 0);
    }
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{

    switch(type){
	case render_type::triangle:
		for(int i = 0; i < state.num_vertices; i+=3){
			data_geometry** tri = new data_geometry*[3];
			for(int j = 0; j < 3; j++){
				tri[j]->data = new float[MAX_FLOATS_PER_VERTEX];
				for(int k = 0; k < state.floats_per_vertex; k++){
					tri[j]->data[k] = state.vertex_data[i+j+k];
				}
			}
			rasterize_triangle(state,(const data_geometry**)tri);
		}
		break;
	case render_type::indexed:
		break;
	case render_type::fan:
		break;
	case render_type::strip:
		break;
	default:
		break;
    }
    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
   data_geometry* out = new data_geometry[3];
   int w = state.image_width;
   int h = state.image_height;
   int i = 0;
   int j = 0;
   unsigned int image_index = 0;

   int ax = 0;
   int ay = 0;
   int bx = 0;
   int by = 0;
   int cx = 0;
   int cy = 0;

   int px = 0;
   int py = 0;

   int AREAabc = 0;
   int AREApbc = 0;
   int AREAapc = 0;
   int AREAabp = 0;

   int alpha = 0;
   int beta = 0;
   int gamma = 0;

   for(int index = 0; index < 3; index++){
	data_vertex* v = 0;
	v->data = in[index]->data;
	const data_vertex* constV = v;
	
	state.vertex_shader(*constV, out[index], state.uniform_data);	

	out[index].gl_Position[0] /= out[index].gl_Position[3];
	out[index].gl_Position[1] /= out[index].gl_Position[3];

	i = w/2.0 * out[index].gl_Position[0] + w/2.0 - (0.5);
	j = h/2.0 * out[index].gl_Position[1] + h/2.0 - (0.5);

	image_index = i + j * w;
	state.image_color[image_index] = make_pixel(255, 255, 255);
   }

   ax = (w/2.0)*out[0].gl_Position[0] + (w/2.0) - (0.5);
   ay = (w/2.0)*out[0].gl_Position[1] + (w/2.0) - (0.5);
   bx = (w/2.0)*out[1].gl_Position[0] + (w/2.0) - (0.5);
   by = (w/2.0)*out[1].gl_Position[1] + (w/2.0) - (0.5);
   cx = (w/2.0)*out[2].gl_Position[0] + (w/2.0) - (0.5);
   cy = (w/2.0)*out[2].gl_Position[1] + (w/2.0) - (0.5);

   AREAabc = 0.5 * (ax * (by - cy)) + (bx * (cy - ay)) + (cx * (ay - by));


   for(px = 0; px < w; px++){
	for(py = 0; py < h; py++){
		AREApbc = 0.5 * (px * (by - cy)) + (bx * (cy - py)) + (cx * (px - by));
   		AREAapc = 0.5 * (ax * (py - cy)) + (px * (cy - ay)) + (cx * (ay - py));
   		AREAabp = 0.5 * (ax * (by - py)) + (bx * (py - ay)) + (px * (ay - by));
		
		alpha = AREApbc / AREAabc;
		beta = AREAapc / AREAabc;
		gamma = AREAabp / AREAabc;

		image_index = px + py * w;

		if(alpha >= 0 && beta >= 0 && gamma >= 0){
			state.image_color[image_index] = make_pixel(255, 255, 255);
		}
	}
   }

   /*AREAabc = 0.5 * (ax * (by - cy)) + (bx * (cy - ay)) + (cx * (ay - by));

   AREApbc = 0.5 * (px * (by - cy)) + (bx * (cy - py)) + (cx * (px - by));
   AREAapc = 0.5 * (ax * (py - cy)) + (px * (cy - ay)) + (cx * (ay - py));
   AREAabp = 0.5 * (ax * (by - py)) + (bx * (py - ay)) + (px * (ay - by));
 
   alpha = AREApbc / AREAabc;
   beta = AREAapc / AREAabc;
   gamma = AREAabp / AREAabc;

   if(alpha >= 0 && beta >= 0 && gamma >= 0){
	state.image_color[image_index] = make_pixel(255, 255, 255);	
   }		
   */
std::cout<<"TODO: implement rasterization"<<std::endl;

}
