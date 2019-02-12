varying vec4 position;
varying vec3 normal;
varying vec3 light_direction;

void main()
{
    /*vec4 ambient = vec4(1, 0, 0, 1);
    vec4 diffuse = vec4(0, 1, 0, 1);
    vec4 specular = vec4(0, 0, 1, 1);*/

    vec4 ldir = vec4(light_direction, 1);
    vec4 fColor = gl_LightModel.ambient;
	
    vec4 norm = normalize(vec4(normal, 0));


    for(int i = 0; i < gl_MaxLights; i++){

    	vec4 light_vecs  = ldir;
    	vec4 view_ray = normalize(-position);
    	vec4 reflection_ray  = -light_vecs + float(2) * dot(normalize(light_vecs), norm) * norm;

    	vec4 ambient = gl_FrontMaterial.ambient * gl_LightSource[i].ambient;
    	vec4 diffuse = gl_FrontMaterial.diffuse * gl_LightSource[i].diffuse * max(dot(norm,normalize(light_vecs)), 0.0);
    	vec4 specular = gl_FrontMaterial.specular * gl_LightSource[i].specular * pow(max(dot(view_ray, reflection_ray),0.0), gl_FrontMaterial.shininess);

	fColor += ambient + diffuse + specular;
    }
    gl_FragColor = fColor;
}
