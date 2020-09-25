#!/usr/bin/env python
# (C) 2018, Tom Eulenfeld, MIT license

from seisgeo import plot_html

def color_func(name):
    return (
        '#3186cc' if name[:2] in ('GP', 'TR') else
        'violet' if name[:2] == 'EL' else
        'green' if name[0] in 'EDOUVW' else
        'red')


plot_html('GEO/geo_eubabrunn_tag3.txt', 'map.html', zone=(33, 'N'), delimiter=',', color_func=color_func)
