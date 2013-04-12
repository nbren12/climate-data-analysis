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


def lat_pressure_contour_cdms(x,tit=None,nlev=None):
    lat = x.getLatitude().getValue()
    lev = x.getLevel().getValue()
    plotme = squeeze(x.getValue())

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

var_names = dict((
    ( 'rlut' , 'toa_net_longwave_flux'),
    ( 'rst'  , 'toa_net_shortwave_flux'),
    ( 'rls'  , 'surface_net_longwave_flux'),
    ( 'rss'  , 'surface_net_shortwave_flux'),
    ( 'hfss' , 'surface_sensible_heat_flux'),
    ( 'hfls' , 'surface_latent_heat_flux'))
)
# Global Energy Budget#{{{
if False:

    var = []
    for key in var_names:
        td = {}
        tt = cdutil.averager(fc02(key),axis="xy",weights="generate")
        cdutil.setTimeBoundsMonthly(tt)
        tt = cdutil.YEAR(tt)
        td = (key,{'val':tt,'name':key,'title':var_names[key]})
        var.append(td)

    var = dict(var)


    # Top of Atmospher Global Energy Budget
    if False:
        toa_lw = var['rlut']['val']
        toa_sw = var['rst']['val']

        net = toa_lw + toa_sw
        fig = figure(figsize=(15,5))
        fig.suptitle('TOA Energy Budget Warm Sun')
        subplot(131)
        plot(toa_lw.getValue())
        title('TOA LW Flux')
        ylabel('W/m^2')
        xlabel('Year')
        subplot(132)
        plot(toa_sw.getValue())
        title('TOA SW Flux')
        xlabel('Year')
        subplot(133)
        plot(net.getValue())
        title('TOA Net Flux')
        xlabel('Year')
        savefig('/Users/noah/Desktop/plot.png',bbox_inches='tight')

    # Surface Energy Budget
    if False:
        rls = var['rls']['val']
        rss = var['rss']['val']
        hfss = var['hfss']['val']
        hfls = var['hfls']['val']

        net = rls + rss + hfss + hfls

        fig = figure(figsize=(8,10))
        fig.suptitle('Surface Budget C02 x 2 ')
        subplot(321)
        plot(rls.copy())
        title('Surface LW')
        subplot(322)
        plot(rss.copy())
        title('Surface SW')
        subplot(323)
        plot(hfss.copy())
        title('Surface Sensible')
        subplot(324)
        plot(hfss.copy())
        title('Surface Latent')
        subplot(325)
        plot (net.copy())
        title('net')

        savefig('/Users/noah/Desktop/plot.png',bbox_inches='tight')

        pass
#}}}

# Latidudinal Energy Budget#{{{
if False:
    var = []
    run_title = 'C02 x 2'
    for key in var_names:
        tt = cdutil.averager(fc02(key),axis="x",weights="generate")
        cdutil.setTimeBoundsMonthly(tt)
        tt = cdutil.YEAR(tt)[-10:]
        tt = cdutil.YEAR.climatology(tt)
        td = (key,tt)
        var.append(td)

    var = dict(var)
    lats = tt.getLatitude()[:]
    gv = lambda tt : squeeze(tt.getValue().copy())
    ave = lambda tt : cdutil.averager(tt,axis="y",weights="generate")

    def Net2Xport(net):
        tmp = net.clone()
        tmp2 = net.clone()
        tmp[tmp<0]=0
        tmp2[tmp2>0]=0
        return min(abs(ave(tmp)),abs(ave(tmp2)))
    # TOA Budget
    rst  = gv(var['rst'])
    rlut = gv(var['rlut'])
    fig = figure(figsize=(6,10))
    subplot(311)
    plot(lats,rst,'r')
    plot(lats,-rlut,'b')

    net = var['rst']+var['rlut']
    net_toa = net.clone()
    msk = gv(net)>0
    fill_between(lats[msk ],rst[msk],-rlut[msk],color=(.8,.8,.8))

    total_trans = Net2Xport(net_toa)
    legend(('Incoming SW','Outgoing LW'),loc=8)
    title('%s TOA Total Xport = %f'%(run_title,total_trans))
    ylabel('Watts per Square Meter')

    # Surface Budget
    rls = var['rls']
    rss = var['rss']
    hfss = var['hfss']
    hfls = var['hfls']

    net_surf = rls + rss + hfss + hfls
    ocean_xport = Net2Xport(net_surf)


    subplot(312)
    plot(lats,gv(net_surf))
    title('%s Net Surface Flux ; Ocean Xport = %.3f'%(run_title,ocean_xport))
    ylabel('Watts per Square Meter')

    # Atmo Budget
    subplot(313)
    atmo_xport = Net2Xport(net_toa-net_surf)
    plot(lats,gv(net_toa-net_surf))
    xlabel('latitude')
    ylabel('Watts per Square Meter')
    title('%s Atmo Budget ; Atmo Xport = %.3f'%(run_title,total_trans-ocean_xport))
    savefig('/Users/noah/Desktop/%s Lat Budget.png'%run_title)
#}}}

if True:

    from mpl_toolkits.basemap import Basemap,addcyclic,shiftgrid,cm

    f = fc02

    net_surf = f('rls')+f('rss')+f('hfss')+f('hfls')
    net_toa  = f('rlut')+f('rst')

    cdutil.setTimeBoundsMonthly(net_surf)
    cdutil.setTimeBoundsMonthly(net_toa)

    def Last10Ave(var):
        return cdutil.YEAR.climatology(cdutil.YEAR(var)[-10:])

    net_surf = Last10Ave(net_surf)
    net_toa = Last10Ave(net_toa)


    var = net_toa
    longitude = var.getLongitude()[:]
    lat = var.getLatitude()[:]
    lon0 = -180

    vargrid,lon= shiftgrid(180.+lon0,squeeze(var.getValue()),longitude,start=False)
    vargrid,lon= addcyclic(vargrid,lon)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.05,0.05,0.9,0.9])
    m = Basemap(projection='robin',lat_0 = 0,lon_0=lon0,resolution='c')
    x, y = m(*np.meshgrid(lon[:], lat[:]))
    im1 = m.pcolor(x,y,vargrid,cmap=plt.cm.RdBu,vmin=-100,vmax=100)
    colorbar()
    m.drawcoastlines()
    # m.fillcontinents()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
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

    zon_p_plots('ta',fc,fc02,title="Double C02 Temperature",nlev=arange(-5,9))
    zon_p_plots('ta',fc,fs0,title="S0 + 2% Temperature",nlev= arange(-5,9))
    zon_p_plots('ua',fc,fc02,title="Double C02 U")
    zon_p_plots('ua',fc,fs0,title="S0 + 2% U")
    zon_p_plots('hus',fc,fc02,title="Double C02 q")
    zon_p_plots('hus',fc,fs0,title="S0 + 2% q")
#}}}
show()
