library(RODBC)
library(rgdal)
library(plotGoogleMaps)
library(qdap)
library(sqldf)
library(igraph)
library(xtable)
library(Hmisc)
library("rattle")
rattle()

clusters_plot <- function(country,centres)
{
  #accepts as inputs AS geolocation data and the number of clusters and calculates
  #clusters and cluster centers
  #outputs a googlmap plot of ASs and cluster centers 
  #outputs dataframe where each AS is assigned a cluster center. 
  #appends cluster centers to the end of the dataframe
  #due to googlemaps limitatation, in a situation where there is an exact coordinates
  #match, only 1 point will be displayed
  set.seed(1)
  y <- cbind(country$max_latitude, country$max_longtitude)
  colnames(y) <- c("latitude", "longtitude")
  cl <- kmeans(y,centres,20,1)
  temp <- cbind(cl$cluster,country)
  centres <- list(colnames(country))
  for(i in seq(1:nrow(cl$center)))
  {
    temp[nrow(temp)+1,10] = cl$center[i,1]
    temp[nrow(temp),12] = cl$center[i,2]
    temp[nrow(temp),11] = cl$center[i,1]
    temp[nrow(temp),13] = cl$center[i,2]
    temp[nrow(temp),1] = i*100
  }
  country_cluster <- temp
  coordinates(temp)<-~min_longtitude+min_latitude
  proj4string(temp)<- CRS("+proj=longlat +datum=WGS84")
  m1<-plotGoogleMaps(temp,filename='myMap1.htm',mapTypeId='ROADMAP',layerName = 'AS cluster and center')
  return (country_cluster)   
}

enrich_country_data <- function(country_connections,country)
{
  #accepts as inputs AS connectivity data and the country in question 
  #looks up clusters of the ASs of the country in question 
  #If not found outputs 'NA'
  #outputs a dataframe with connectivity data enriched with clusters connectivity
  country_connections <- sqlQuery(channel,country_connections)
  country_connections_with_cluster <- sqldf("SELECT country_connections.*,country.cl_cluster as cluster2 from country_connections left join country on country.id=country_connections.id2")
  country_all <- sqldf("SELECT country.cl_cluster as cluster1,country_connections_with_cluster.* from country_connections_with_cluster left join country on country.id=country_connections_with_cluster.id")
  country_final <- sqldf("SELECT cluster1,id,reg_country,name,IP_prefix_counter,city,country,id2,relationship,reg_country_1 as reg_country2,name_1 as name2,cluster2 from country_all")
  return (country_final)
}

create_cluster_graph <- function (country_final,number_of_clusters)
{
  #accepts as inputs country final data file and number of clusters and outputs
  #a graph of country cluster connectivity
  total_clusters = number_of_clusters + 1   
  country_graph <- graph.empty() + vertices(1:total_clusters)
  for(i in seq(1:nrow(country_final)))
  {
    #print(i)
    from <- country_final[i,1]
    to <- country_final[i,12]
    if(is.na(as.character(to))) 
      to = number_of_clusters+1
    if(is.na(country_final[i,9])==FALSE)
    {
      if(country_final[i,9] == 'c2p')
      {
        country_graph <- add.edges(country_graph,c(from,to))
      }
      if(country_final[i,9] == 'p2p')
      {
        country_graph <- add.edges(country_graph,c(from,to))
        country_graph <- add.edges(country_graph,c(to,from))
      }
    }  
  }
  
  return (country_graph)
}

calculate_metrics <- function (country_cluster_graph, final_frame)
{
  #Accepts country cluster graph and country name and outputs 
  #degree, betweenness and closeness metrics 
  metrics <- matrix()
  degree <- degree(country_cluster_graph)
  betweenness <- betweenness(country_cluster_graph)
  closeness <- closeness(country_cluster_graph)
  number_of_clusters<-max(final_frame['cluster1']) 
  number_of_as_in_cluster<-numeric(number_of_clusters) 
  inbound <- numeric(number_of_clusters)
  outbound <- numeric(number_of_clusters)
  for (i in seq(1:number_of_clusters)) 
  { 
    number_of_as_in_cluster[i] <- nrow(unique(final_frame[final_frame['cluster1']==i,]['id']))
    outbound[i] <- fn$sqldf("Select count(*) from final_frame b where b.cluster1 = $i and (b.relationship = 'c2p' or b.relationship = 'p2p')")
    inbound[i] <- fn$sqldf("Select count(*) from final_frame b where b.cluster1 = $i and (b.relationship = 'p2c' or b.relationship = 'p2p')")
  }
    number_of_as_in_cluster[number_of_clusters+1] <- NA
    outbound[number_of_clusters+1] <- NA
    inbound[number_of_clusters+1] <- NA  
  metrics <- cbind(degree,betweenness,closeness,outbound,inbound)    
  return (metrics)
}

belarus_cluster_metrics <- calculate_metrics(belarus_cluster_graph,belarus_final)


calculate_as_per_cluster<-function(final_frame,cluster_col,id_col) 
{ 
  number_of_clusters<-max(final_frame[cluster_col]) 
  vec<-numeric(number_of_clusters) 
  for (i in seq(1:number_of_clusters)) 
  { vec[i]<-nrow(unique(final_frame[final_frame[cluster_col]==i,][id_col])) } 
  return(vec) 
}

descriptive_stats <- function(country_all)
  #accepts all connections within a country and calculates basic descriptive statistics
{
  stats <- list()
#  stats['COUNTRY CODE'] <- sqldf("Select distinct reg_country from country_all")
  stats['ASs'] <- sqldf("Select count(*) from (Select distinct b.id from country_all b) as t1")
  stats['Connections'] <- sqldf("Select count(*) from country_all c")
  stats['BIG AS'] <- sqldf("Select count(*) from country_all c where c.distance <> 0")  
  stats['SMALL AS'] <- sqldf("Select count(*) from country_all c where c.distance = 0")  
  stats['BIG AS ABROAD'] <- sqldf("Select count(*) from country_all c where c.distance <> 0 and c.reg_country_1 <> c.reg_country")  
  stats['SMALL AS ABROAD'] <- sqldf("Select count(*) from country_all c where c.distance = 0 and c.reg_country_1 <> c.reg_country")
  stats['BIG AS HOME'] <- sqldf("Select count(*) from country_all c where c.distance <> 0 and c.reg_country_1 = c.reg_country")  
  stats['SMALL AS HOME'] <- sqldf("Select count(*) from country_all c where c.distance = 0 and c.reg_country_1 = c.reg_country")
  return (stats)
}

#CONSTANTS
BELARUS_NUMBER_OF_CLUSTERS = 3
RUSSIA_NUMBER_OF_CLUSTERS = 13
CHINA_NUMBER_OF_CLUSTERS = 6
GERMANY_NUMBER_OF_CLUSTERS = 6
FRANCE_NUMBER_OF_CLUSTERS = 9
JAPAN_NUMBER_OF_CLUSTERS = 6
SAUDI_NUMBER_OF_CLUSTERS = 2
UKRAINE_NUMBER_OF_CLUSTERS = 8
US_NUMBER_OF_CLUSTERS = 26
GB_NUMBER_OF_CLUSTERS = 6
ES_NUMBER_OF_CLUSTERS = 7
IT_NUMBER_OF_CLUSTERS = 6
CA_NUMBER_OF_CLUSTERS = 6
BR_NUMBER_OF_CLUSTERS = 10
KR_NUMBER_OF_CLUSTERS = 5
ZA_NUMBER_OF_CLUSTERS = 4
AR_NUMBER_OF_CLUSTERS = 5

#Database channel definition
channel <- odbcConnect("AS", uid="root", pwd="c1edd62dad")

#SQL selects choosing ASs that are both registered and geo-located in one country only
China_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'CN' and g.country = 'CN') GROUP BY a.id HAVING distance = 0;"
Germany_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'DE' and g.country = 'DE') GROUP BY a.id HAVING distance = 0;"
France_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'FR' and g.country = 'FR') GROUP BY a.id HAVING distance = 0;"
Russia_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'RU' and g.country = 'RU') GROUP BY a.id HAVING distance = 0;"
Ukraine_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'UA' and g.country = 'UA') GROUP BY a.id HAVING distance = 0;"
Japan_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'JP' and g.country = 'JP') GROUP BY a.id HAVING distance = 0;"
Belarus_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'BY' and g.country = 'BY') GROUP BY a.id HAVING distance = 0;"
Saudi_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'SA' and g.country = 'SA') GROUP BY a.id HAVING distance = 0;"
US_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'US' and g.country = 'US') GROUP BY a.id HAVING distance = 0;"

#SQL selects choosing all connections of ASs in a particular country irrespective of the distance within them
China_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'CN' and g.country = 'CN')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Germany_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'DE' and g.country = 'DE')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
France_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'FR' and g.country = 'FR')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Russia_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'RU' and g.country = 'RU')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Ukraine_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'UA' and g.country = 'UA')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Japan_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'JP' and g.country = 'JP')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Belarus_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'BY' and g.country = 'BY')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Saudi_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'SA' and g.country = 'SA')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
US_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'US' and g.country = 'US')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"


#SQL lookups of the connections (id connected to, type of relationship) of the ASs in one country only
China_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'CN' and g.country = 'CN')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Germany_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'DE' and g.country = 'DE')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
France_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'FR' and g.country = 'FR')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Russia_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'RU' and g.country = 'RU')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Ukraine_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'UA' and g.country = 'UA')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Japan_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'JP' and g.country = 'JP')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Belarus_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'BY' and g.country = 'BY')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
Saudi_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'SA' and g.country = 'SA')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"
US_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'US' and g.country = 'US')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

#Data regarding AS geo-location retrieval
China_one <- sqlQuery(channel,China_one_location)
Germany_one <- sqlQuery(channel,Germany_one_location)
France_one <- sqlQuery(channel,France_one_location)
Russia_one <- sqlQuery(channel,Russia_one_location)
Ukraine_one <- sqlQuery(channel,Ukraine_one_location)
Japan_one <- sqlQuery(channel,Japan_one_location)
Belarus_one <- sqlQuery(channel,Belarus_one_location)
Saudi_one <- sqlQuery(channel,Saudi_one_location)
US_one <- sqlQuery(channel,US_one_location)

#Data regarding AS geo-location retrieval
China_all <- sqlQuery(channel,China_all_connections)
Germany_all <- sqlQuery(channel,Germany_all_connections)
France_all <- sqlQuery(channel,France_all_connections)
Russia_all <- sqlQuery(channel,Russia_all_connections)
Ukraine_all <- sqlQuery(channel,Ukraine_all_connections)
Japan_all <- sqlQuery(channel,Japan_all_connections)
Belarus_all <- sqlQuery(channel,Belarus_all_connections)
Saudi_all <- sqlQuery(channel,Saudi_all_connections)
US_all <- sqlQuery(channel,US_all_connections)

#Calculating clusters and outputing maps
belarus <- clusters_plot(Belarus_one,BELARUS_NUMBER_OF_CLUSTERS)
russia <- clusters_plot(Russia_one,RUSSIA_NUMBER_OF_CLUSTERS)
china <- clusters_plot(China_one,CHINA_NUMBER_OF_CLUSTERS)
germany <- clusters_plot(Germany_one,GERMANY_NUMBER_OF_CLUSTERS)
france <- clusters_plot(France_one,FRANCE_NUMBER_OF_CLUSTERS)
japan <- clusters_plot(Japan_one,JAPAN_NUMBER_OF_CLUSTERS)
saudi <- clusters_plot(Saudi_one,SAUDI_NUMBER_OF_CLUSTERS)
ukraine <- clusters_plot(Ukraine_one,UKRAINE_NUMBER_OF_CLUSTERS)
us <- clusters_plot(US_one,US_NUMBER_OF_CLUSTERS)

#Enriching data with connectivity and clusters information
belarus_final <- enrich_country_data(Belarus_connected_to,belarus)
russia_final <- enrich_country_data(Russia_connected_to,russia)
ukraine_final <- enrich_country_data(Ukraine_connected_to,ukraine)
china_final <- enrich_country_data(China_connected_to,china)
germany_final <- enrich_country_data(Germany_connected_to,germany)
france_final <- enrich_country_data(France_connected_to,france)
japan_final <- enrich_country_data(Japan_connected_to,japan)
saudi_final <- enrich_country_data(Saudi_connected_to,saudi)
us_final <- enrich_country_data(US_connected_to,us)


#Creating country cluster graphs
belarus_cluster_graph <- create_cluster_graph(belarus_final,BELARUS_NUMBER_OF_CLUSTERS)
russia_cluster_graph <- create_cluster_graph(russia_final,RUSSIA_NUMBER_OF_CLUSTERS)
ukraine_cluster_graph <- create_cluster_graph(ukraine_final,UKRAINE_NUMBER_OF_CLUSTERS)
china_cluster_graph <- create_cluster_graph(china_final,CHINA_NUMBER_OF_CLUSTERS)
germany_cluster_graph <- create_cluster_graph(germany_final,GERMANY_NUMBER_OF_CLUSTERS)
france_cluster_graph <- create_cluster_graph(france_final,FRANCE_NUMBER_OF_CLUSTERS)
japan_cluster_graph <- create_cluster_graph(japan_final,JAPAN_NUMBER_OF_CLUSTERS)
saudi_cluster_graph <- create_cluster_graph(saudi_final,SAUDI_NUMBER_OF_CLUSTERS)
us_cluster_graph <- create_cluster_graph(us_final,US_NUMBER_OF_CLUSTERS)

#calculating metrics for country cluster graphs
belarus_cluster_metrics <- calculate_metrics(belarus_cluster_graph,belarus_final)
russia_cluster_metrics <- calculate_metrics(russia_cluster_graph,russia_final)
ukraine_cluster_metrics <- calculate_metrics(ukraine_cluster_graph,ukraine_final)
china_cluster_metrics <- calculate_metrics(china_cluster_graph,china_final)
germany_cluster_metrics <- calculate_metrics(germany_cluster_graph,germany_final)
france_cluster_metrics <- calculate_metrics(france_cluster_graph,france_final)
japan_cluster_metrics <- calculate_metrics(japan_cluster_graph,japan_final)
saudi_cluster_metrics <- calculate_metrics(saudi_cluster_graph,saudi_final)
us_cluster_metrics <- calculate_metrics(us_cluster_graph,us_final)

#Country stats
Belarus_stats <- descriptive_stats(Belarus_all)
Russia_stats <- descriptive_stats(Russia_all)
Ukraine_stats <- descriptive_stats(Ukraine_all)
China_stats <- descriptive_stats(China_all)
Germany_stats <- descriptive_stats(Germany_all)
France_stats <- descriptive_stats(France_all)
Japan_stats <- descriptive_stats(Japan_all)
Saudi_stats <- descriptive_stats(Saudi_all)
US_stats <- descriptive_stats(US_all)

stats_merged <- rbind(Belarus_stats,Russia_stats,Ukraine_stats,China_stats,Germany_stats,France_stats,Japan_stats,Saudi_stats,US_stats)

#AS with most IP prefixes in a country
Belarus_data_by_prefix <- sqldf("Select * from Belarus_all b group by b.id order by IP_prefix_counter DESC")
Russia_data_by_prefix <- sqldf("Select * from Russia_all b group by b.id order by IP_prefix_counter DESC")
Ukraine_data_by_prefix <- sqldf("Select * from Ukraine_all b group by b.id order by IP_prefix_counter DESC")
China_data_by_prefix <- sqldf("Select * from China_all b group by b.id order by IP_prefix_counter DESC")
Germany_data_by_prefix <- sqldf("Select * from Germany_all b group by b.id order by IP_prefix_counter DESC")
France_data_by_prefix <- sqldf("Select * from France_all b group by b.id order by IP_prefix_counter DESC")
Japan_data_by_prefix <- sqldf("Select * from Japan_all b group by b.id order by IP_prefix_counter DESC")
Saudi_data_by_prefix <- sqldf("Select * from Saudi_all b group by b.id order by IP_prefix_counter DESC")
US_data_by_prefix <- sqldf("Select * from US_all b group by b.id order by IP_prefix_counter DESC")

#GB
GB_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'GB' and g.country = 'GB') GROUP BY a.id HAVING distance = 0;"

GB_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'GB' and g.country = 'GB')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

GB_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'GB' and g.country = 'GB')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

GB_one <- sqlQuery(channel,GB_one_location)

GB_all <- sqlQuery(channel,GB_all_connections)

GB <- clusters_plot(GB_one,GB_NUMBER_OF_CLUSTERS)

GB_final <- enrich_country_data(GB_connected_to,GB)

GB_cluster_graph <- create_cluster_graph(GB_final,GB_NUMBER_OF_CLUSTERS)

GB_cluster_metrics <- calculate_metrics(GB_cluster_graph,GB_final)

GB_stats <- descriptive_stats(GB_all)

#SPAIN

ES_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'ES' and g.country = 'ES') GROUP BY a.id HAVING distance = 0;"

ES_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'ES' and g.country = 'ES')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

ES_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'ES' and g.country = 'ES')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

ES_one <- sqlQuery(channel,ES_one_location)

ES_all <- sqlQuery(channel,ES_all_connections)

ES <- clusters_plot(ES_one,ES_NUMBER_OF_CLUSTERS)

ES_final <- enrich_country_data(ES_connected_to,ES)

ES_cluster_graph <- create_cluster_graph(ES_final,ES_NUMBER_OF_CLUSTERS)

ES_cluster_metrics <- calculate_metrics(ES_cluster_graph,ES_final)

ES_stats <- descriptive_stats(ES_all)

#ITALY

IT_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'IT' and g.country = 'IT') GROUP BY a.id HAVING distance = 0;"

IT_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'IT' and g.country = 'IT')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

IT_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'IT' and g.country = 'IT')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

IT_one <- sqlQuery(channel,IT_one_location)

IT_all <- sqlQuery(channel,IT_all_connections)

IT <- clusters_plot(IT_one,IT_NUMBER_OF_CLUSTERS)

IT_final <- enrich_country_data(IT_connected_to,IT)

IT_cluster_graph <- create_cluster_graph(IT_final,IT_NUMBER_OF_CLUSTERS)

IT_cluster_metrics <- calculate_metrics(IT_cluster_graph,IT_final)

IT_stats <- descriptive_stats(IT_all)

#CANADA

CA_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'CA' and g.country = 'CA') GROUP BY a.id HAVING distance = 0;"

CA_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'CA' and g.country = 'CA')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

CA_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'CA' and g.country = 'CA')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

CA_one <- sqlQuery(channel,CA_one_location)

CA_all <- sqlQuery(channel,CA_all_connections)

CA <- clusters_plot(CA_one,CA_NUMBER_OF_CLUSTERS)

CA_final <- enrich_country_data(CA_connected_to,CA)

CA_cluster_graph <- create_cluster_graph(CA_final,CA_NUMBER_OF_CLUSTERS)

CA_cluster_metrics <- calculate_metrics(CA_cluster_graph,CA_final)

CA_stats <- descriptive_stats(CA_all)

#BRAZIL

BR_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'BR' and g.country = 'BR') GROUP BY a.id HAVING distance = 0;"

BR_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'BR' and g.country = 'BR')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

BR_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'BR' and g.country = 'BR')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

BR_one <- sqlQuery(channel,BR_one_location)

BR_all <- sqlQuery(channel,BR_all_connections)

BR <- clusters_plot(BR_one,BR_NUMBER_OF_CLUSTERS)

BR_final <- enrich_country_data(BR_connected_to,BR)

BR_cluster_graph <- create_cluster_graph(BR_final,BR_NUMBER_OF_CLUSTERS)

BR_cluster_metrics <- calculate_metrics(BR_cluster_graph,BR_final)

BR_stats <- descriptive_stats(BR_all)

#SOUTH KOREA

KR_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'KR' and g.country = 'KR') GROUP BY a.id HAVING distance = 0;"

KR_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'KR' and g.country = 'KR')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

KR_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'KR' and g.country = 'KR')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

KR_one <- sqlQuery(channel,KR_one_location)

KR_all <- sqlQuery(channel,KR_all_connections)

KR <- clusters_plot(KR_one,KR_NUMBER_OF_CLUSTERS)

KR_final <- enrich_country_data(KR_connected_to,KR)

KR_cluster_graph <- create_cluster_graph(KR_final,KR_NUMBER_OF_CLUSTERS)

KR_cluster_metrics <- calculate_metrics(KR_cluster_graph,KR_final)

KR_stats <- descriptive_stats(KR_all)

#ARGENTINA

AR_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'AR' and g.country = 'AR') GROUP BY a.id HAVING distance = 0;"

AR_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'AR' and g.country = 'AR')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

AR_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'AR' and g.country = 'AR')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

AR_one <- sqlQuery(channel,AR_one_location)

AR_all <- sqlQuery(channel,AR_all_connections)

AR <- clusters_plot(AR_one,AR_NUMBER_OF_CLUSTERS)

AR_final <- enrich_country_data(AR_connected_to,AR)

AR_cluster_graph <- create_cluster_graph(AR_final,AR_NUMBER_OF_CLUSTERS)

AR_cluster_metrics <- calculate_metrics(AR_cluster_graph,AR_final)

AR_stats <- descriptive_stats(AR_all)

#SOUTH AFRICA

ZA_one_location <- "SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country, max(g.latitude) as max_latitude, min(g.latitude) as min_latitude, max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude, sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2)) as distance FROM IP_GEO as g LEFT JOIN IP_AS as i ON g.ip=i.ip JOIN AS.AS as a ON i.as_id=a.id WHERE (a.reg_country = 'ZA' and g.country = 'ZA') GROUP BY a.id HAVING distance = 0;"

ZA_all_connections <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'ZA' and g.country = 'ZA')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

ZA_connected_to <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'ZA' and g.country = 'ZA')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

ZA_one <- sqlQuery(channel,ZA_one_location)

ZA_all <- sqlQuery(channel,ZA_all_connections)

ZA <- clusters_plot(ZA_one,ZA_NUMBER_OF_CLUSTERS)

ZA_final <- enrich_country_data(ZA_connected_to,ZA)

ZA_cluster_graph <- create_cluster_graph(ZA_final,ZA_NUMBER_OF_CLUSTERS)

ZA_cluster_metrics <- calculate_metrics(ZA_cluster_graph,ZA_final)

ZA_stats <- descriptive_stats(ZA_all)

#UKRAINE 2014

Ukraine_all_connections_2014 <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'UA' and g.country = 'UA')  GROUP BY a.id) as t1 LEFT JOIN AS_connectivity_2014 as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

Ukraine_connected_to_2014 <- "SELECT * FROM (SELECT a.*,count(a.id) as IP_prefix_counter, g.city, g.country,  max(g.latitude) as max_latitude, min(g.latitude) as min_latitude,  max(g.longtitude) as max_longtitude, min(g.longtitude) as min_longtitude,  sqrt(pow(sqrt(pow(max(g.latitude)-min(g.latitude),2)),2)  +pow(sqrt(pow(max(g.longtitude)-min(g.longtitude),2)),2))  as distance  FROM IP_GEO as g  LEFT JOIN IP_AS as i  ON g.ip=i.ip  JOIN AS.AS as a  ON i.as_id=a.id  WHERE (a.reg_country = 'UA' and g.country = 'UA')  GROUP BY a.id  HAVING distance = 0) as t1 LEFT JOIN AS_connectivity_2014 as c ON t1.id = c.id LEFT JOIN AS.AS as a on c.id2 = a.id;"

Ukraine_all_2014 <- sqlQuery(channel,Ukraine_all_connections_2014)

Ukraine_stats_2014 <- descriptive_stats(Ukraine_all_2014)

Ukraine_data_by_prefix_2014 <- sqldf("Select * from Ukraine_all_2014 b group by b.id order by IP_prefix_counter DESC")

ukraine_final_2014 <- enrich_country_data(Ukraine_connected_to_2014,ukraine)



IP_AS_GEO <- sqlQuery(channel,"SELECT * FROM IP_AS JOIN IP_GEO ON IP_AS.IP = IP_GEO.IP")


