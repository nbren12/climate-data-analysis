from pylab import *
import netCDF4 as nc
from ipdb import set_trace as st

# Functions#{{{
def monthAve(x):
    s = list( x.shape )
    s[0] = 12
    out = zeros(s)
    for i in range(12):
        out[i,:,:] = x[i::12,:,:].mean(axis=0)
    return out

def seasonAve(x):
    pass
def lat_pressure_contour(lat,lev,plotme,tit=None,nlev = None):
    if nlev is None:
        nlev = 9
    fig = figure(figsize = (8,6))
    contourf(lat,lev,plotme,nlev)
    gca().invert_yaxis()
    colorbar()
    cs = contour(lat,lev,plotme,nlev,colors='k')
    clabel(cs,inline=1)
    xlabel('Latitude')
    ylabel('Pressure (hPa)')
    xticks(arange(-80,90,20))
    axis('tight')

    if tit is not None:
        title(tit)

    return fig
#}}}

# Load Data
fc02 = nc.Dataset('dc02.nc')
fs0  = nc.Dataset('warmsun.nc')
fc  = nc.Dataset('control.nc')
fstd = nc.Dataset('standard.nc')


lat = fc02.variables['lat'][:]
lon = fc02.variables['lon'][:]
lev = fc02.variables['lev'][:]
lats = lat.shape[0]
lons = lon.shape[0]

# Global Energy Budget
from mpl_toolkits.basemap import Basemap
# create figure, axes instances.
var_names = dict((
    ( 'rlut' , 'toa_net_longwave_flux'),
    ( 'rst'  , 'toa_net_shortwave_flux'),
    ( 'rls'  , 'surface_net_longwave_flux'),
    ( 'rss'  , 'surface_net_shortwave_flux'),
    ( 'hfss' , 'surface_sensible_heat_flux'),
    ( 'hfls' , 'surface_latent_heat_flux'))
)
cV = fc.variables['rst']
cVz = cV[:].mean(axis=-1)

import cdutil


fig = plt.figure(figsize=(8,6))
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='robin',lat_0 = 90,lon_0=lon.mean(),resolution='c')
x, y = m(*np.meshgrid(lon, lat))
im1 = m.pcolor(x,y,fc02.variables['zg'][0,2,:,:],cmap=plt.cm.jet)
m.drawcoastlines()
# m.fillcontinents()
m.drawmapboundary()
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
# Response Diagnostic#{{{
if False:
    def zon_p_plots(var,dsetA,dsetB,title="",nlev=None):

        controlV = dsetA.variables[var]
        compareV = dsetB.variables[var]

        # Average Zonally over last 20 years
        bt = 30*12

        aV = mean(controlV[bt:,:,:,:],axis=(-1))
        bV = mean(compareV[bt:,:,:,:],axis=(-1))


        aVm = monthAve(aV)
        bVm = monthAve(bV)

        cVm = bVm - aVm

        plotme = cVm[[11,0,1],:,:].mean(axis=0)
        lat_pressure_contour(lat,lev,plotme,tit='%s DJF'%title,nlev=nlev)

        plotme = cVm[6:9,:,:].mean(axis=0)
        lat_pressure_contour(lat,lev,plotme,tit='%s JJA'%title,nlev=nlev)

    zon_p_plots('ta',fc,fc02,title="Double C02 Temperature",nlev=arange(-5,9))
    zon_p_plots('ta',fc,fs0,title="S0 + 2% Temperature",nlev= arange(-5,9))
    zon_p_plots('ua',fc,fc02,title="Double C02 U")
    zon_p_plots('ua',fc,fs0,title="S0 + 2% U")
    zon_p_plots('hus',fc,fc02,title="Double C02 q")
    zon_p_plots('hus',fc,fs0,title="S0 + 2% q")
#}}}

show()

