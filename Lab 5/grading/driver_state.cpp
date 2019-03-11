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
    state.image_depth = new float[length];
    for(int i = 0; i < length; i++){
	state.image_depth[i] = 1;
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
	
	// call the vertex shader on each vertex
	//
	
/*	data_geometry dg[num_vertices];
	float dg_data[num_vertices*floats_per_vertex];

	for(  each vertex)
{	
	data_vertex v;
	v.data = vertex_data[ right spot];

	dg[i].data = dg_data [ right spot];

	vertex_shader(); 
}*/

    switch(type){
	case render_type::triangle: {
		const data_geometry* dgarr[3];
		data_geometry dg[3];
		data_vertex dv[3];
		int num_tri = state.num_vertices / 3;
		int index = 0;
		
		for(int i = 0; i < num_tri; i++){
			for(int j = 0; j < 3; j++){
				dv[j].data = &state.vertex_data[index];
				dg[j].data = dv[j].data;
				state.vertex_shader(dv[j], dg[j], state.uniform_data);
				dgarr[j] = &dg[j];
				index += state.floats_per_vertex;
			}
		//	rasterize_triangle(state, dgarr);
			clip_triangle(state, dgarr, 0);
		}
		break;
	}		
	case render_type::indexed:{
		const data_geometry* dgarr[3];
		data_geometry dg[3];
		data_vertex dv[3];

		for(int i = 0; i < 3 * state.num_triangles; i += 3){
			for(int j = 0; j < 3; j++){
				dv[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
				dg[j].data = dv[j].data;
				state.vertex_shader(dv[j], dg[j], state.uniform_data);
				dgarr[j] = &dg[j];
			}
			clip_triangle(state, dgarr, 0);
		}
		break;
	}
	case render_type::fan: {
		const data_geometry* dgarr[3];
		data_geometry dg[3];
		data_vertex dv[3];

		for(int i = 0; i < state.num_vertices; i++){
			for(int j = 0; j < 3; j++){
				int index = i + j;
				if(j == 0){
					index = 0;
				}
				dv[j].data = &state.vertex_data[index * state.floats_per_vertex];
				dg[j].data = dv[j].data;
				state.vertex_shader(dv[j], dg[j], state.uniform_data);
				dgarr[j] = &dg[j];
			}
			clip_triangle(state, dgarr, 0);
		}
		break;
	}
	case render_type::strip: {
		const data_geometry* dgarr[3];
		data_geometry dg[3];
		data_vertex dv[3];

		for(int i = 0; i < state.num_vertices - 2; i++){
			for(int j = 0; j < 3; j++){
				dv[j].data = &state.vertex_data[(i + j) * state.floats_per_vertex];
				dg[j].data = dv[j].data;
				state.vertex_shader(dv[j], dg[j], state.uniform_data);
				dgarr[j] = &dg[j];
			}
			clip_triangle(state, dgarr, 0);
		}
		break;
	}
	default: {
		break;
	}
    }
  //  std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==1)
    {
        rasterize_triangle(state, in);
        return;
    }
  //  std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    
    vec4 a = in[0]->gl_Position;
    vec4 b = in[1]->gl_Position;
    vec4 c = in[2]->gl_Position;

    const data_geometry *input[3] = {in[0], in[1], in[2]};
    data_geometry new_data0[3];
    data_geometry new_data1[3];
    float a0, a1, b0, b1;
    vec4 p0, p1;
 
    if (a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3])
        return;
    else 
        if (a[2] < -a[3] && b[2] >= -b[3] && c[2] >= -c[3])
        {
            b0 = (-b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
            b1 = (-a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);
            p0 = b0 * a + (1 - b0) * b;
            p1 = b1 * c + (1 - b1) * a;

            new_data0[0].data = new float[state.floats_per_vertex];
            new_data0[1] = *in[1];
            new_data0[2] = *in[2];

            for (int i = 0; i < state.floats_per_vertex; ++i)
                switch (state.interp_rules[i])
                {
                case interp_type::flat:
                    new_data0[0].data[i] = in[0]->data[i];
                    break;
                case interp_type::smooth:
                    new_data0[0].data[i] = b1 * in[2]->data[i] + (1 - b1) * in[0]->data[i];
                    break;
                case interp_type::noperspective:
                    a0 = b1 * in[2]->gl_Position[3] / (b1 * in[2]->gl_Position[3] + (1 - b1) * in[0]->gl_Position[3]);
                    new_data0[0].data[i] = a0 * in[2]->data[i] + (1 - a0) * in[0]->data[i];
                    break;
                default:
                    break;
                }

            new_data0[0].gl_Position = p1;
            input[0] = &new_data0[0];
            input[1] = &new_data0[1];
            input[2] = &new_data0[2];
            
            clip_triangle(state, input, face + 1);

            new_data1[0].data = new float[state.floats_per_vertex];
            new_data1[2] = *in[2];

            for (int i = 0; i < state.floats_per_vertex; ++i)
                switch (state.interp_rules[i])
                {
                case interp_type::flat:
                    new_data1[0].data[i] = in[0]->data[i];
                    break;
                case interp_type::smooth:
                    new_data1[0].data[i] = b0 * in[0]->data[i] + (1 - b0) * in[1]->data[i];
                    break;
                case interp_type::noperspective:
                    a0 = b0 * in[0]->gl_Position[3] / (b0 * in[0]->gl_Position[3] + (1 - b0) * in[1]->gl_Position[3]);
                    new_data1[0].data[i] = a0 * in[0]->data[i] + (1 - a0) * in[1]->data[i];
                    break;
                default:
                    break;
                }

            new_data1[0].gl_Position = p0;
            input[0] = &new_data1[0];
            input[1] = &new_data0[1];
            input[2] = &new_data0[0];
        }

    clip_triangle(state, input, face + 1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
   int w = state.image_width;
   int h = state.image_height;

   float ax = 0;
   float ay = 0;
   float bx = 0;
   float by = 0;
   float cx = 0;
   float cy = 0;

   int px = 0;
   int py = 0;

   float AREAabc = 0;
   float AREApbc = 0;
   float AREAapc = 0;
   float AREAabp = 0;

   float alpha = 0;
   float beta = 0;
   float gamma = 0;
   
   float alphp = 0;
   float betp = 0;
   float gamp = 0;

   ax = (w/2.0)*(in[0]->gl_Position[0] / in[0]->gl_Position[3]) + (w/2.0) - (0.5);
   ay = (h/2.0)*(in[0]->gl_Position[1] / in[0]->gl_Position[3]) + (h/2.0) - (0.5);
   bx = (w/2.0)*(in[1]->gl_Position[0] / in[1]->gl_Position[3]) + (w/2.0) - (0.5);
   by = (h/2.0)*(in[1]->gl_Position[1] / in[1]->gl_Position[3]) + (h/2.0) - (0.5);
   cx = (w/2.0)*(in[2]->gl_Position[0] / in[2]->gl_Position[3]) + (w/2.0) - (0.5);
   cy = (h/2.0)*(in[2]->gl_Position[1] / in[2]->gl_Position[3]) + (h/2.0) - (0.5);

   float minX = std::min(ax, std::min(bx, cx));
   float minY = std::min(ay, std::min(by, cy));
   float maxX = std::max(ax, std::max(bx, cx));
   float maxY = std::max(ay, std::max(by, cy));

   //std::cout << ax << " " << ay << " " << bx << " " << by << " " << cx << " " << cy << std::endl;

    AREAabc = 0.5 * (((bx*cy)-(cx*by)) - ((ax*cy)-(cx*ay)) + ((ax*by)-(bx*ay)));
   // AREAabc = (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
   //std::cout << AREAabc << std::endl;

   if(minX < 0){
	minX = 0;
   }
   if(maxX > state.image_width){
	maxX = state.image_width;
   }
   if(minY < 0){
	minY = 0;
   }
   if(maxY > state.image_height){
	maxY = state.image_height;
   }

   for(px = minX; px < maxX; px++){
	for(py = minY; py < maxY; py++){
		int index = px + py * state.image_width;
		AREApbc = 0.5 * (((bx*cy) - (cx*by)) + ((by-cy)*px) + ((cx-bx)*py));
   		//AREApbc = (px - bx) * (cy - by) - (py - by) * (cx - bx);
		AREAapc = 0.5 * (((cx*ay) - (ax*cy)) + ((cy-ay)*px) + ((ax-cx)*py));
		//AREAapc = (px - cx) * (ay - cy) - (py - cy) * (ax - cx);
		AREAabp = 0.5 * (((ax*by) - (bx*ay)) + ((ay-by)*px) + ((bx-ax)*py));
		//AREAabp = (px - ax) * (by - ay) - (py - ay) * (bx - ax);		

		alpha = AREApbc / AREAabc;
		beta = AREAapc / AREAabc;
		gamma = AREAabp / AREAabc;


//		std::cout << alpha << " " << beta << " " << gamma << std::endl;		

		if(alpha >= 0 && beta >= 0 && gamma >= 0){
			data_fragment frag;
			frag.data = new float[MAX_FLOATS_PER_VERTEX];
			data_output o;
			float z_dep = (alpha * (in[0]->gl_Position[2] / in[0]->gl_Position[3])) + (beta * (in[1]->gl_Position[2] / in[1]->gl_Position[3])) + (gamma * (in[2]->gl_Position[2] / in[2]->gl_Position[3]));
			if(state.image_depth[index] > z_dep){
				for(int findex = 0; findex < state.floats_per_vertex; findex++){
					float fl;
					switch(state.interp_rules[findex]){
						case interp_type::flat:
							frag.data[findex] = in[0]->data[findex];
							break;
						case interp_type::smooth:
							fl = (alpha/in[0]->gl_Position[3]) +(beta/in[1]->gl_Position[3]) + (gamma/in[2]->gl_Position[3]);
							alphp = alpha /  (in[0]->gl_Position[3] * fl);
							betp = beta / (in[1]->gl_Position[3] * fl);
							gamp = gamma / (in[2]->gl_Position[3] * fl);
							frag.data[findex] = (alphp * in[0]->data[findex]) + (betp * in[1]->data[findex]) + (gamp * in[2]->data[findex]);
							break;
						case interp_type::noperspective:
							frag.data[findex] = alpha*in[0]->data[findex] + beta*in[1]->data[findex] + gamma*in[2]->data[findex];
							break;
						default:
							break;
					}
				}
				state.fragment_shader(frag, o, state.uniform_data);
			
				state.image_color[index] = make_pixel(o.output_color[0] * 255, o.output_color[1] * 255, o.output_color[2] * 255);
				state.image_depth[index] = z_dep;
			}
		}
	}
//std::cout<<"TODO: implement rasterization"<<std::endl;
    }
}
