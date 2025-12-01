clc;
clear all;
#load trabalho_missao_.mat;



# Coordenadas
lat = - [15 52 16 ]; # s
long = - [47 55 07]; # w
h = [1066];

geodesica_ini = [lat, long, h];

h_dec = h;
lat_dec = lat(1)+lat(2)/(60)+lat(3)/(60*60);
long_dec = long(1)+long(2)/(60)+long(3)/(60*60);

# Conversao Geodesica -> Cartesiano
a = 6378137;
e = 0.0818;
N0 = a/(sqrt(1-e^2*(sind(lat_dec))^2));
rx(1) = (N0+h_dec)*cosd(lat_dec)*cosd(long_dec);
ry(1) = (N0+h_dec)*cosd(lat_dec)*sind(long_dec);
rz(1) = (N0*(1-e^2)+h_dec)*sind(lat_dec);

#Determina??o das velocidades iniciais
vx(1) = 0;
vy(1) = 0;
vz(1) = 0;

# Determina??o da atitude inicial
phi(1) = atitude_ini(1);
theta(1) = atitude_ini(2);
psi(1) = atitude_ini(3);

# Descobrir a rota??o da atitude
m=length(ax_ms2);

for i = 1 : m
  #Determina??o da matriz de rota??o do gir?metro
  Rb2n = [ 1 sind(phi(i))*tand(theta(i)) cosd(phi(i))*tand(theta(i));
           0        cosd(phi(i))                -sind(phi(i));
           0 sind(phi(i))/cosd(theta(i)) cosd(phi(i))/cosd(theta(i))];

  #Determina??o das derivadas dos ?ngulos de atitude
  d_ang = Rb2n * [wx_degs(i);
                  wy_degs(i);
                  wz_degs(i)];

  #Determina??o dos ?ngulos phi e theta segundo os gir?metros
  phi(i+1) = phi(i) + d_ang(1) * dt;
  theta(i+1) = theta(i) + d_ang(2) * dt;

  #Determina??o do psi com o mag
  Cb2n = [cosd(theta(i)) sind(phi(i))*sind(theta(i)) cosd(phi(i))*sind(theta(i));
          0                       cosd(phi(i))             -sind(phi(i));
          -sind(theta(i)) sind(phi(i))*cosd(theta(i)) cosd(phi(i))*cosd(theta(i))]; #matriz de rota??o para os gir?metros com psi=0

  #Determina??o dos ms para encontrar o psi
  m = Cb2n * [mx_uT(i);
              my_uT(i);
              mz_uT(i)];

  #Determina??o do psi verdadeiro com a declina??o
  psi(i+1) = atan2d(-m(2) , m(1)) + declina;

  #Determina??o das atitudes finais com os accs
  Rn2b = [cosd(psi(i))*cosd(theta(i)) cosd(psi(i))*sind(theta(i))*sind(phi(i))-sind(psi(i))*cosd(phi(i)) cosd(psi(i))*sind(theta(i))*cosd(phi(i))+sind(psi(i))*sind(phi(i));
          sind(psi(i))*cosd(theta(i)) sind(psi(i))*sind(theta(i))*sind(phi(i))+cosd(psi(i))*cosd(phi(i)) sind(psi(i))*sind(theta(i))*cosd(phi(i))-cosd(psi(i))*sind(phi(i));
          -sind(theta(i))             cosd(theta(i))*sind(phi(i))                                        cosd(theta(i))*cosd(phi(i))]; #matriz de rota??o para os aceler?metros

  #Determina??o das acelera??es
  a = Rn2b * [ax_ms2(i);
              ay_ms2(i);
              az_ms2(i)];

  #Retirando a gravidade do az
  a(3) = a(3) - g(i);

  #Determina??o das velocidades
  vx(i+1) = vx(i) + a(1) * dt;
  vy(i+1) = vy(i) + a(2) * dt;
  vz(i+1) = vz(i) + a(3) * dt;

  #Determina??o das posi??es
  rx(i+1) = rx(i) + vx(i+1) * dt;
  ry(i+1) = ry(i) + vy(i+1) * dt;
  rz(i+1) = rz(i) + vz(i+1) * dt;
  endfor




