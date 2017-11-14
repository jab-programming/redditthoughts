import os
import sys
import time
import math
import itertools
import xml.dom.minidom

from nmeaSentence import makeAltobridgeSentence

class ScenarioFeeder:
    def __init__(self,scenario):
        self.scenario = scenario
        self.suffix = ','.join([self.scenario.origin,self.scenario.destination])
        self.start_time = time.time()

    def getLine(self):
        seconds_since_start = time.time() - self.start_time
        latitude, longitude, altitude = self.scenario.position_after(seconds_since_start)
        nmea = makeAltobridgeSentence(latitude,longitude,'%4.4f' % altitude)
        return '%s;%s' % ( nmea, self.suffix )

    def atEnd(self):
        return False

def knots_to_metres_per_second(knots):
    '''1 knot = 0.000514444444 km/s

    According to http://www.google.ie/#q=1+knots+in+m%2Fs
    '''
    return knots * 0.514444444

class ScenarioLeg:
    def __init__(self, start, end ):
        self.start = start
        self.end = end
        self.metres_per_second = knots_to_metres_per_second(start[3])
        self.metres = 1000 * vincenty_kilometres(start[0], start[1], end[0], end[1])
        self.seconds = int( self.metres / self.metres_per_second )

    def __str__(self):
        if self.metres < 1000:
            return "%8.2f m from (%6.2f,%7.2f) at %6.2f to (%6.2f,%7.2f) at %6.2f" % ( self.metres, self.start[0], self.start[1], self.start[2], self.end[0], self.end[1], self.end[2] )
        else:
            return "%8.2fkm from (%6.2f,%7.2f) at %6.2f to (%6.2f,%7.2f) at %6.2f" % ( self.metres / 1000.0, self.start[0], self.start[1], self.start[2], self.end[0], self.end[1], self.end[2] )

    def __repr__(self):
        return '<ScenarioLeg "%s">' % self

    def stages(self,seconds_per_point):
        times = self.seconds / seconds_per_point
        assert float(times) == int(times)
        def points(start, end):
            diff = end - start
            step = diff / times
            result = [ ]
            next = start
            for i in range(0,times):
                result += [ next ]
                next += step
            return result
        latitudes = points(self.start[0], self.end[0])
        longitudes = points(self.start[1], self.end[1])
        altitudes = points(self.start[2], self.end[2])
        return [ x for x in  itertools.izip(latitudes, longitudes, altitudes) ]

    def position_after(self,seconds_after_start):
        ratio = float(seconds_after_start) / self.seconds
        assert 0.0 <= ratio <= 1.0
        # print self.seconds, seconds_after_start, ratio
        start_latitude, start_longitude, start_altitude, start_knots = self.start
        end_latitude, end_longitude, end_altitude, end_knots = self.end
        latitude = start_latitude + ((end_latitude - start_latitude) * ratio)
        longitude = start_longitude + ((end_longitude - start_longitude) * ratio)
        altitude = start_altitude + ((end_altitude - start_altitude) * ratio)
        return latitude, longitude, altitude

class Scenario:
    def add_turns(self, turns):
        self.turns = turns
        self.legs = []
        previous = turns[0]
        first = previous
        for turn in self.turns[1:]:
            self.legs += [ ScenarioLeg(previous, turn) ]
            previous = turn
        self.legs += [ ScenarioLeg(turn, first) ]

    def __str__(self):
        return '%d turns from %s to %s' % ( len(self.turns), self.origin, self.destination )

    def stages(self,seconds_per_point):
        result = []
        for leg in self.legs:
            result += leg.stages(seconds_per_point)
        return result

    def position_after(self,seconds_after_start):
        remainder = seconds_after_start
        i = 0
        while True:
            for leg in self.legs:
                seconds_after_leg_start = remainder
                remainder -= leg.seconds
                if remainder < 0:
        #            print 'Leg', i
                    return leg.position_after(seconds_after_leg_start)
                i += 1
        return -999, -999


equatorial_radius_in_metres = 6378.137
polar_radius_in_metres = 6356.7523142
acceptable_accuracy_in_metres =  0.00001

def law_of_cosines_kilometres(start_latitude, start_longitude, end_latitude, end_longitude ):
    start_latitude = math.radians(start_latitude)
    end_latitude = math.radians(end_latitude)
    latitude_sin = math.sin(start_latitude) * math.sin(end_latitude)
    latitude_cos = math.cos(start_latitude) * math.cos(end_latitude)

    longitude_difference = math.radians(end_longitude - start_longitude)
    difference_cos = math.cos(longitude_difference)

    average_latitude = (start_latitude + end_latitude) / 2.0
    earthRadius = equatorial_radius_in_metres - (21 * math.sin(average_latitude))
    return earthRadius * math.acos(latitude_sin + latitude_cos * difference_cos)

def haversine_kilometres(start_latitude, start_longitude, end_latitude, end_longitude):
    latitude_difference = math.radians(end_latitude - start_latitude)
    longitude_difference = math.radians(end_longitude - start_longitude)
    latitude_h = math.sin(latitude_difference/2.0) ** 2.0
    longitude_h = math.sin(longitude_difference/2.0) ** 2.0
    start_latitude = math.radians(start_latitude)
    end_latitude = math.radians(end_latitude)
    a = latitude_h + math.cos(start_latitude) * math.cos(end_latitude) * longitude_h
    return equatorial_radius_in_metres * 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))

def vincenty_kilometres(start_latitude, start_longitude, end_latitude, end_longitude):
    start_latitude = math.radians(start_latitude)
    end_latitude = math.radians(end_latitude)
    start_longitude = math.radians(start_longitude)
    end_longitude = math.radians(end_longitude)

    f = 1 / 298.257223563
    longitudeDifference = end_longitude - start_longitude
    U1 = math.atan((1 - f) * math.tan(start_latitude))
    U2 = math.atan((1 - f) * math.tan(end_latitude))
    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    fred = longitudeDifference
    previous_fred = 2 * math.pi
    iterLimit = 20
    cosSqAlpha = 0.0
    sinSigma = 0.0
    cosSigma = 0.0
    cos2SigmaM = 0.0
    sigma = 0.0
    while math.fabs(fred - previous_fred) > acceptable_accuracy_in_metres and iterLimit > 0:
        sin_fred = math.sin(fred)
        cos_fred = math.cos(fred)
        sinSigma = math.sqrt((cosU2 * sin_fred) * (cosU2 * sin_fred) + (cosU1 * sinU2 - sinU1 * cosU2 * cos_fred) * (cosU1 * sinU2 - sinU1 * cosU2 * cos_fred))
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cos_fred
        sigma = math.atan2(sinSigma, cosSigma)
        alpha = math.asin(cosU1 * cosU2 * sin_fred / sinSigma)
        cosSqAlpha = math.cos(alpha) * math.cos(alpha)
        cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        previous_fred = fred
        fred = longitudeDifference + (1 - C) * f * math.sin(alpha) * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)))
        iterLimit -= 1
    if not iterLimit:
        raise ValueError("Vincenty failed to converge")
    radiusToEquatorSquared = equatorial_radius_in_metres * equatorial_radius_in_metres
    radiusToPoleSquared = polar_radius_in_metres * polar_radius_in_metres
    uSq = cosSqAlpha * (radiusToEquatorSquared - radiusToPoleSquared) / (radiusToPoleSquared)
    radiusToEquator = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    radiusToPole = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = radiusToPole * sinSigma * (cos2SigmaM + radiusToPole / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - radiusToPole / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)))
    return polar_radius_in_metres * radiusToEquator * (sigma - deltaSigma)

def read_scenario_xml(path_to_xml_file,node_name):
    if not os.path.isfile(path_to_xml_file):
        return []
    dom = xml.dom.minidom.parse(path_to_xml_file)
    top_node = dom.getElementsByTagName(node_name).item(0)
    return top_node.childNodes

def parse_scenarios(nodes):
    values = {}
    for node in nodes:
        if node.hasChildNodes():
            value = node.childNodes[0].nodeValue
            values[node.nodeName] = values.get(node.nodeName,[]) + [ value.strip() ]
    for key, value in values.iteritems():
        if type(value) == type([]) and len(value) == 1:
            values[key] = value[0]
    longitudes = []
    latitudes = []
    if 'arcViewData' in values:
        lines = values['arcViewData'].splitlines()
        assert lines[0] == 'Polyline'
        for line in lines:
            try: counter, longitude, latitude, _, _ = line.split()
            except ValueError: continue
            longitudes += [ float(longitude) ]
            latitudes += [ float(latitude) ]
        assert lines[-1] == 'END'
    if 'altitude' in values:
        altitudes = [ float(a) for a in values['altitude'] ]
    if 'speed' in values:
        knots = [ float(s) for s in values['speed'] ]
    scenario = Scenario()
    scenario.add_turns( [ i for i in itertools.izip(latitudes, longitudes, altitudes, knots) ] )
    scenario.origin = values['origin']
    scenario.destination = values['destination']
    scenario.time_ratio = values['timeRatio']
    return scenario


def read_scenario_file(path_to_scenario_file):
    nodes = read_scenario_xml(path_to_scenario_file,'scenarioSetup');
    if not nodes:
        return None
    return parse_scenarios(nodes)

def run_scenario(path_to_scenario):
    scenario = read_scenario_file(path_to_scenario)
    feeder = ScenarioFeeder(scenario)
    while not feeder.atEnd():
        print feeder.getLine()
        time.sleep(1)

def main(args):
    for arg in args:
        if os.path.isfile(arg):
            run_scenario(arg)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
