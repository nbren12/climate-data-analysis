from pylab import *
import netCDF4 as nc
from ipdb import set_trace as st
import cdms2 as cdms
import cdutil

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

def lat_pressure_contour_cdms(x,tit=None,nlev=None):
    lat = x.getLatitude().getValue()
    lev = x.getLevel().getValue()
    plotme = squeeze(x.getValue())
    st()

    lat_pressure_contour(lat,lev,plotme,tit=tit,nlev=nlev)

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

fc02 = cdms.open('dc02.nc')
fs0  = cdms.open('warmsun.nc')
fc  = cdms.open('control.nc')
fstd = cdms.open('standard.nc')

# Global Energy Budget
var_names = dict((
    ( 'rlut' , 'toa_net_longwave_flux'),
    ( 'rst'  , 'toa_net_shortwave_flux'),
    ( 'rls'  , 'surface_net_longwave_flux'),
    ( 'rss'  , 'surface_net_shortwave_flux'),
    ( 'hfss' , 'surface_sensible_heat_flux'),
    ( 'hfls' , 'surface_latent_heat_flux'))
)

var = []
for key in var_names:
    td = {}
    tt = cdutil.averager(fc02(key),axis="xy",weights="generate")
    td = (key,{'val':tt,'name':key,'title':var_names[key]})
    var.append(td)

var = dict(var)
figure(figsize=(8,6))
net = var['rlut']['val']+var['rst']['val']
plot(net.getValue())

# lat = fc02.variables['lat'][:]
# lon = fc02.variables['lon'][:]
# lev = fc02.variables['lev'][:]
# lats = lat.shape[0]
# lons = lon.shape[0]

# Global Energy Budget
# create figure, axes instances.



# from mpl_toolkits.basemap import Basemap
# fig = plt.figure(figsize=(8,6))
# ax = fig.add_axes([0.05,0.05,0.9,0.9])
# m = Basemap(projection='robin',lat_0 = 90,lon_0=lon.mean(),resolution='c')
# x, y = m(*np.meshgrid(lon, lat))
# im1 = m.pcolor(x,y,fc02.variables['zg'][0,2,:,:],cmap=plt.cm.jet)
# m.drawcoastlines()
# # m.fillcontinents()
# m.drawmapboundary()
# m.drawparallels(np.arange(-90.,120.,30.))
# m.drawmeridians(np.arange(0.,420.,60.))
# cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
# Response Diagnostic#{{{
if False:
    def zon_p_plots(var,dsetA,dsetB,title="",nlev=None):
        # Average Zonally over last 20 years
        bt = 30*12

        controlV = dsetA(var,time=slice(bt,600))
        compareV = dsetB(var,time=slice(bt,600))

        a = cdutil.averager(controlV,axis="x")
        b = cdutil.averager(compareV,axis="x")

        c =b-a
        cdutil.setTimeBoundsMonthly(c)
        plotme = cdutil.DJF.climatology(c)
        lat_pressure_contour_cdms(plotme,tit='%s DJF'%title,nlev=nlev)

        plotme = cdutil.JJA.climatology(c)
        lat_pressure_contour_cdms(plotme,tit='%s JJA'%title,nlev=nlev)

    c = zon_p_plots('ta',fc,fc02,title="Double C02 Temperature",nlev=arange(-5,9))
    zon_p_plots('ta',fc,fs0,title="S0 + 2% Temperature",nlev= arange(-5,9))
    zon_p_plots('ua',fc,fc02,title="Double C02 U")
    zon_p_plots('ua',fc,fs0,title="S0 + 2% U")
    zon_p_plots('hus',fc,fc02,title="Double C02 q")
    zon_p_plots('hus',fc,fs0,title="S0 + 2% q")
#}}}

show()

