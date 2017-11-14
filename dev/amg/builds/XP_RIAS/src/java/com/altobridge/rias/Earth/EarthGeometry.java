/**
 * Copyright (c) 2010 Altobridge, Ireland. All Rights Reserved
 * THIS SOURCE CODE IS CONFIDENTIAL AND PROPRIETARY AND MAY NOT
 * BE USED OR DISTRIBUTED WITHOUT THE WRITTEN PERMISSION OF
 * Altobridge.
 *
 */
package com.altobridge.rias.Earth;

import com.altobridge.rias.util.RiasException;

import java.util.logging.Logger;

public class EarthGeometry {

    private static final double EQUATORIAL_RADIUS = 6378.137;   // wgs-84 equatorial radius - used in distance calculation
    private static final double POLAR_RADIUS = 6356.7523142;
    private static final double ACCEPTABLE_ACCURACY = 0.001; // value in km, i.e. 1 metres

    protected static Logger logger = Logger.getLogger(EarthGeometry.class.getPackage().getName());
    static private int haversineCount = 0;
    static private int vincentyCount = 0;
    static public void reset()
    {
        haversineCount = 0;
        vincentyCount = 0;

    }
    /**
     * Calculate the distance based on two sets of longitude and latitude coordinates.
     * The coordinates must firstly be converted from degrees to radians for the algorithm to work.
     * The algorithm is "acos[sin(latitude A) * sin(latitude B) + cos(latitude A) * cos(latitude B) * cos(+longitude B - (+longitude A))] * Radius of Earth in kilometres."
     *
     * @param aLatitudeDegrees  - latitude of point A
     * @param aLongitudeDegrees - longitute of point A
     * @param bLatitudeDegrees  - latitude of point B
     * @param bLongitudeDegrees - longitude of point B
     *
     * @return the calculated distance
     *
     */
    private static double haversineDistance(double aLatitudeDegrees, double aLongitudeDegrees, double bLatitudeDegrees, double bLongitudeDegrees)
    {
        ++haversineCount;
        double aLatitude = Math.toRadians(aLatitudeDegrees);
        double bLatitude = Math.toRadians(bLatitudeDegrees);
        final double latitudeSinProduct = Math.sin(aLatitude) * Math.sin(bLatitude);
        final double latitudeCosProduct = Math.cos(aLatitude) * Math.cos(bLatitude);

        double differenceInLongitudeDegrees = bLongitudeDegrees - aLongitudeDegrees;
        double longitudeDifference = Math.toRadians(differenceInLongitudeDegrees);
        final double cosOfDifference = Math.cos(longitudeDifference);

        final double averageLatitude = (aLatitude + bLatitude) / 2.0;
        double earthRadius = EQUATORIAL_RADIUS - (21 * Math.sin(averageLatitude));
        return earthRadius * Math.acos(latitudeSinProduct + latitudeCosProduct * cosOfDifference);
    }

    public static void log(final String prefix)
    {
        logger.info(prefix + " Haversine: " + haversineCount + ", Vincenty: " + vincentyCount);
    }
    public static double greatCircleDistance(double aLatitudeDegrees, double aLongitudeDegrees, double bLatitudeDegrees, double bLongitudeDegrees)
    {
        try {
            return vincentyDistance(aLatitudeDegrees, aLongitudeDegrees, bLatitudeDegrees, bLongitudeDegrees);
        } catch ( ConvergenceException e ) {
            return haversineDistance(aLatitudeDegrees, aLongitudeDegrees, bLatitudeDegrees, bLongitudeDegrees);
        }
    }

    public static double greatCircleDistance(MapPosition a, MapPosition b)
    {
        return greatCircleDistance(a.getLatitude(), a.getLongitude(), b.getLatitude(), b.getLongitude());
    }


    /**
     * A more accurate calculation of Great Circle Distance, which allows better for actual geometry of the earth.
     * It is more computationally expensive.
     * But experience from testing suggests that (given the ACCEPTABLE_ACCURACY of 0.1m) the main loop only iterates twice.
     *
     * @param aLatitudeDegrees  Latitude of point A
     * @param aLongitudeDegrees Longitude of point A
     * @param bLatitudeDegrees  Latitude of point B
     * @param bLongitudeDegrees Longitude of point
     *
     * @return the Great Circle Distance between the two points.
     *
     * @throws ConvergenceException if the loop fails to converge.
     */
    private static double vincentyDistance(double aLatitudeDegrees, double aLongitudeDegrees, double bLatitudeDegrees, double bLongitudeDegrees) throws ConvergenceException
    {
        ++vincentyCount;
        // http://www.movable-type.co.uk/scripts/LatLongVincenty.html

        if ( aLatitudeDegrees == bLatitudeDegrees && aLongitudeDegrees == bLongitudeDegrees )
            return 0.0;

        double aLatitude = Math.toRadians(aLatitudeDegrees);
        double aLongitude = Math.toRadians(aLongitudeDegrees);
        double bLatitude = Math.toRadians(bLatitudeDegrees);
        double bLongitude = Math.toRadians(bLongitudeDegrees);

        double f = 1 / 298.257223563;
        double longitudeDifference = bLongitude - aLongitude;
        double U1 = Math.atan((1 - f) * Math.tan(aLatitude));
        double U2 = Math.atan((1 - f) * Math.tan(bLatitude));
        double sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);
        double sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);
        double lambda = longitudeDifference, previousLambda = 2 * Math.PI;
        int iterLimit = 20;
        double cosSqAlpha = 0.0;
        double sinSigma = 0.0;
        double cosSigma = 0.0;
        double cos2SigmaM = 0.0;
        double sigma = 0.0;
        while ( Math.abs(lambda - previousLambda) > ACCEPTABLE_ACCURACY && --iterLimit > 0 ) {
            double sinLambda = Math.sin(lambda), cosLambda = Math.cos(lambda);
            sinSigma = Math.sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
            sigma = Math.atan2(sinSigma, cosSigma);
            double alpha = Math.asin(cosU1 * cosU2 * sinLambda / sinSigma);
            cosSqAlpha = Math.cos(alpha) * Math.cos(alpha);
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
            double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
            previousLambda = lambda;
            lambda = longitudeDifference + (1 - C) * f * Math.sin(alpha) * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
        }
        if ( iterLimit == 0 ) {
            throw new ConvergenceException("Vincenty failed to converge");
        }
        double radiusToEquatorSquared = EQUATORIAL_RADIUS * EQUATORIAL_RADIUS;
        double radiusToPoleSquared = POLAR_RADIUS * POLAR_RADIUS;
        double uSq = cosSqAlpha * (radiusToEquatorSquared - radiusToPoleSquared) / (radiusToPoleSquared);
        double radiusToEquator = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
        double radiusToPole = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
        double deltaSigma = radiusToPole * sinSigma * (cos2SigmaM + radiusToPole / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - radiusToPole / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        return POLAR_RADIUS * radiusToEquator * (sigma - deltaSigma);

    }

    /**
     * Convert the given direction to a value within (PI,-PI)
     *
     * @param directionInDegrees An angle (from equator) in degrees.
     *
     * @return The angle in radians.
     */
    public static double normalisedRadians(double directionInDegrees)
    {
        double degrees = normalisedDegrees(directionInDegrees);
        return Math.toRadians(degrees);
    }

    public static double normalisedDegrees(double directionInDegrees)
    {
        double degrees = directionInDegrees % 360.0;
        if ( degrees > 180 )
            degrees -= 360.0;
        else if ( degrees < -180 )
            degrees += 360;
        else if ( degrees == -180 )
            return 180;
        return degrees;
    }

    /**
     * Compute the dot product ab . bc
     */
    private static double dotProduct(MapPosition a, MapPosition b, MapPosition c)
    {
        double ab_latitude = b.getLatitude() - a.getLatitude();
        double ab_longitude = b.getLongitude() - a.getLongitude();
        double bc_latitude = c.getLatitude() - b.getLatitude();
        double bc_longitude = c.getLongitude() - b.getLongitude();
        return ab_longitude * bc_longitude + ab_latitude * bc_latitude;
    }

    /**
     * Calculate the distance from a point to a line segment.
     * If the point falls outside the line segment then the
     *     value returned is the distance to the closest end point,
     *     otherwise the "perpendicular" distance from the point to the line is returned.
     *
     * @param start The start point of the line
     * @param end The end point of the line
     * @param point The point whose distance from the line is sought
     */
    public static double distanceLineSegmentToPoint(MapPosition start, MapPosition end, MapPosition point)
    {
        double dot = dotProduct(start,end,point);
        if ( dot > 0 )
            return greatCircleDistance(end,point);
        dot = dotProduct(end,start,point);
        if ( dot > 0 )
            return greatCircleDistance(start,point);
        return distanceLineToPoint(start,end,point);
    }


    /**
     * Calculate the earth-spherical distance from a point to a line.
     * The distance is analogous to the "perpendicular distance" from point to line on a plane,
     *     but allows for Earth curvature.
     */
    private static double distanceLineToPoint(MapPosition start, MapPosition end, MapPosition point)
    {
        // http://mathforum.org/library/drmath/view/51785.html

        MapUnitVector startVector = new MapUnitVector(start);
        MapUnitVector endVector = new MapUnitVector(end);
        MapUnitVector pointVector = new MapUnitVector(point);
        MapUnitVector crossVector = endVector.crossProduct(startVector);
        crossVector.normalise();

        double dot = pointVector.dotProduct(crossVector);
        double angle =  Math.acos(dot);
        angle = Math.PI / 2.0 - angle;
        final double averageLatitude = (start.getLatitude() + end.getLatitude() + point.getLatitude()) / 3.0;
        double earthRadius = EQUATORIAL_RADIUS - (21 * Math.sin(averageLatitude));
        double distance = angle * earthRadius;
        double absoluteDistance = Math.abs(distance);
        return absoluteDistance;
    }
}
