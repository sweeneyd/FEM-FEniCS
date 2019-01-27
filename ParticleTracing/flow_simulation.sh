echo "[!] Generating mesh..."
simulation_file=$1
gmsh_file=$2
mesh_file="${gmsh_file/.geo/.msh}"
xml_file="${gmsh_file/.geo/.xml}"
velocity_file="velocity.h5"
echo "[+] Converting $gmsh_file to $previous_file"
/Applications/Gmsh.app/Contents/MacOS/Gmsh "$gmsh_file" -0 || break
./fenics.sh "dolfin-convert $mesh_file $xml_file || break"
./fenics.sh "python3 $simulation_file $xml_file || break"
./fenics.sh "python3 particle_flow.py $xml_file $velocity_file 25 9 0.1"
