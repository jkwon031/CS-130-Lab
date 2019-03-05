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
				index += state.floats_per_vertex;
			}
			rasterize_triangle(state, dgarr);
		}
//		for(int i = 0; i < state.num_vertices; i += 3){
//	data_geometry* tri[3];
//	tri[0] = &dg[i+0];
			
//data_geometry** tri = new data_geometry*[3];
	//		data_geometry* tmp = new data_geometry[3];
	//		data_geometry** triangle = new data_geometry*[3];

		//	tri[0] = new data_geometry;
		//	tri[1] = new data_geometry;
		//	tri[2] = new data_geometry;


		//	const_cast<data_geometry*>(tri[0])->data = new float[MAX_FLOATS_PER_VERTEX];
		//	const_cast<data_geometry*>(tri[1])->data = new float[MAX_FLOATS_PER_VERTEX];
		//	const_cast<data_geometry*>(tri[2])->data = new float[MAX_FLOATS_PER_VERTEX];
			
		/*	for(int j = 0; j < 3; j++){

				tri[j] = new data_geometry;
				data_vertex v;
				v.data = new float[MAX_FLOATS_PER_VERTEX];
				tri[j]->data = new float[MAX_FLOATS_PER_VERTEX];
				for(int k = 0; k < state.floats_per_vertex; k++){
					tri[0]->data[j] = state.vertex_data[j + (state.floats_per_vertex * i)];
					tri[1]->data[j] = state.vertex_data[j + (state.floats_per_vertex * (i + 1))];
					tri[2]->data[j] = state.vertex_data[j + (state.floats_per_vertex * (i + 2))];*/

					//for(int k = 0; k < 3; k++){
					//const_cast<data_geometry*>(tmp[0])->data = new float[MAX_FLOATS_PER_VERTEX];
					//const_cast<data_geometry*>(tmp[1])->data = new float[MAX_FLOATS_PER_VERTEX];
					//const_cast<data_geometry*>(tmp[2])->data = new float[MAX_FLOATS_PER_VERTEX];
					
			//		v.data[k] = state.vertex_data[k + (state.floats_per_vertex * (i + j))];
					//tri[j]->data[k] = v.data[k];
			//	}
			//	state.vertex_shader((const data_vertex)v, *tri[j], state.uniform_data);
			//}

			//state.vertex_shader(tri[0], tmp[0], state.uniform_data);
			//state.vertex_shader(tri[1], tmp[1], state.uniform_data);
			//state.vertex_shader(tri[2], tmp[2], state.uniform_data);

			//for(int k = 0; k < 3; k++){
				//triangle[i] = &tmp[i];
			//}

		//	rasterize_triangle(state,(const data_geometry**)tri);
	//	}
		/*	delete [] tri[0]->data;
			delete [] tri[1]->data;
			delete [] tri[2]->data;
			delete tri[0];
			delete tri[1];
			delete tri[2];
			delete [] tri;
			delete [] tmp;
			//delete [] triangle[0]->data;
			//delete [] triangle[1]->data;
			//delete [] triangle[2]->data;
			delete triangle[0];
			delete triangle[1];
			delete triangle[2];
			delete [] triangle;*/
		break;
	}		
	case render_type::indexed:{
		break;
	}
	case render_type::fan: {
		break;
	}
	case render_type::strip: {
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
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
  //  std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
 //  data_geometry* out;
   float w = state.image_width;
   float h = state.image_height;

//   int i = 0;
  // int j = 0;

   float ax = 0;
   float ay = 0;
   float bx = 0;
   float by = 0;
   float cx = 0;
   float cy = 0;

   float px = 0;
   float py = 0;

   float AREAabc = 0;
   float AREApbc = 0;
   float AREAapc = 0;
   float AREAabp = 0;

   float alpha = 0;
   float beta = 0;
   float gamma = 0;
   
 //  for(int index = 0; index < 3; index++){
	//out[index] = in[index];

	//out[index].gl_Position[0] = in[index]->gl_Position[0] / in[index]->gl_Position[3];
	//out[index].gl_Position[1] = in[index]->gl_Position[1] / in[index]->gl_Position[3];

//	i = w/2.0 * out[index].gl_Position[0] + w/2.0 - (0.5);
//	j = h/2.0 * out[index].gl_Position[1] + h/2.0 - (0.5);

//	int image_index = i + j * w;
	//state.image_color[image_index] = make_pixel(255, 255, 255);
 //  }



 /*vec3 w = {in[0]->gl_Position[3], in[1]->gl_Position[3], in[2]->gl_Position[3]};
    vec3 x = {in[0]->gl_Position[0] / w[0], in[1]->gl_Position[0] / w[1], in[2]->gl_Position[0] / w[2]};
    vec3 y = {in[0]->gl_Position[1] / w[0], in[1]->gl_Position[1] / w[1], in[2]->gl_Position[1] / w[2]};
    vec3 z = {in[0]->gl_Position[2] / w[0], in[1]->gl_Position[2] / w[1], in[2]->gl_Position[2] / w[2]};

    vec2 v0 = {(w / 2) * x[0] + (w / 2) - 0.5, (h / 2) * y[0] + (h / 2) - 0.5};
    vec2 v1 = {(float)((state.image_width / 2) * x[1]) + (float)((state.image_width / 2) - 0.5), 
               (float)((state.image_height / 2) * y[1]) + (float)((state.image_height / 2) - 0.5)};
    vec2 v2 = {(float)((state.image_width / 2) * x[2]) + (float)((state.image_width / 2) - 0.5), 
               (float)((state.image_height / 2) * y[2]) + (float)((state.image_height / 2) - 0.5)};

    vec3 area = {(v2[0] - v0[0]) * (v1[1] - v0[1]) - (v2[1] - v0[1]) * (v1[0] - v0[0]), 0, 0};
*/


   ax = (w/2.0)*(in[0]->gl_Position[0] / in[0]->gl_Position[3]) + (w/2.0) - (0.5);
   ay = (h/2.0)*(in[0]->gl_Position[1] / in[0]->gl_Position[3]) + (h/2.0) - (0.5);
   bx = (w/2.0)*(in[1]->gl_Position[0] / in[1]->gl_Position[3]) + (w/2.0) - (0.5);
   by = (h/2.0)*(in[1]->gl_Position[1] / in[1]->gl_Position[3]) + (h/2.0) - (0.5);
   cx = (w/2.0)*(in[2]->gl_Position[0] / in[2]->gl_Position[3]) + (w/2.0) - (0.5);
   cy = (h/2.0)*(in[2]->gl_Position[1] / in[2]->gl_Position[3]) + (h/2.0) - (0.5);

   //std::cout << ax << " " << ay << " " << bx << " " << by << " " << cx << " " << cy << std::endl;

   //AREAabc = 0.5 * (ax * (by - cy)) + (bx * (cy - ay)) + (cx * (ay - by));
   AREAabc = 0.5 * (((bx*cy) - (cx*by))-((ax*cy) - (cx*ay)) - ((ax*by)-(bx*ay)));

   //std::cout << AREAabc << std::endl;


   for(px = 0; px < w; px++){
	for(py = 0; py < h; py++){
		//AREApbc = 0.5 * (((bx*cy) - (cx*by)) - ((px*cy) - (cx*py)) - ((px*by)-(bx*py)));
		AREApbc = 0.5 * (((bx*cy) - (cx*by)) + ((by-cy)*px) + ((cx-bx)*py));
   		//AREAapc = 0.5 * (((px*cy) - (cx*py)) - ((ax*cy) - (cx*ay)) - ((ax*py)-(px*ay)));
		AREAapc = 0.5 * (((cx*ay) - (ax*cy)) + ((cy-ay)*px) + ((ax-cx)*py));
   		//AREAabp = 0.5 * (((bx*py) - (px*by)) - ((ax*py) - (px*ay)) - ((ax*by)-(bx*ay)));
		AREAabp = 0.5 * (((ax*by) - (bx*ay)) + ((ay-by)*px) + ((bx-ax)*py));		

		alpha = AREApbc / AREAabc;
		beta = AREAapc / AREAabc;
		gamma = AREAabp / AREAabc;


		std::cout << alpha << " " << beta << " " << gamma << std::endl;		

		if(alpha >= 0 && beta >= 0 && gamma >= 0){
			//state.image_color[image_index] = make_pixel(255, 255, 255);
			int image_index = px + py * w;
			auto *data = new float[MAX_FLOATS_PER_VERTEX];
			data_fragment frag{data};
			data_output o;

			for(int findex = 0; findex < state.floats_per_vertex; findex++){
				float fl;
				switch(state.interp_rules[findex]){
					case interp_type::flat:
						frag.data[findex] = in[0]->data[findex];
						break;
					case interp_type::smooth:
						fl = (alpha/in[0]->gl_Position[3] + beta/in[1]->gl_Position[3] + gamma/in[2]->gl_Position[3]);
						alpha /= (fl * (in[0]->gl_Position[3]));
						beta /= (fl * (in[1]->gl_Position[3]));
						gamma /= (fl * (in[2]->gl_Position[3]));
						break;
					case interp_type::noperspective:
						frag.data[findex] = alpha*in[0]->data[findex] + beta*in[1]->data[findex] + gamma*in[2]->data[findex];
						break;
					default:
						break;
				}
			}
			state.fragment_shader((const data_fragment)frag, o, state.uniform_data);
			o.output_color *= 255;

			state.image_color[image_index] = make_pixel(o.output_color[0], o.output_color[1], o.output_color[2]);
		}
	}
   }
//std::cout<<"TODO: implement rasterization"<<std::endl;

}
