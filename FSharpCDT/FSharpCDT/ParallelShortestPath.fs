open System
open System.Threading.Tasks
open System.Collections.Generic
open System.Collections.Concurrent

type City =
    | Boise     | LosAngeles    | NewYork   | Seattle
    | StLouis   | Phoenix       | Boston    | Chicago
    | Denver

// Known distances

let distTweenCities =
    let distances = new Dictionary<City, Dictionary<City, int>>()
    [
        (Boise, Seattle, 496);      (Boise,   Denver,  830);    
        (Boise,   Chicago, 1702);   (Seattle, LosAngeles, 1141); 
        (Seattle, Denver,  1321);   (LosAngeles, Denver,  1022);  
        (LosAngeles, Phoenix, 371); (Phoenix, Denver,  809);      
        (Phoenix, StLouis, 1504);   (Denver,  StLouis, 8588);     
        (Denver,  Chicago, 1009);   (Chicago, NewYork, 811);     
        (Chicago, Boston,  986);    (StLouis, Chicago, 300);
        (Boston, StLouis,  986);    (NewYork, Boston,  211)
    ]

    // Convert list of links between cities into a dictionary
    |> List.iter (fun (cityA, cityB, dist) ->

        if not <| distances.ContainsKey(cityA) then
            distances.Add(cityA, new Dictionary<City, int>())
        if not <| distances.ContainsKey(cityB) then
            distances.Add(cityB, new Dictionary<City, int>())

        distances.[cityA].Add(cityB, dist)
        distances.[cityB].Add(cityA, dist))

     // Return our dictionary of dictionary of distances
     distances

let stdout = Console.Out

let shortestPathBetween startingCity finalDestination =

    // Keep track of the shortest path from startCity to a given city and
    // the path used to get to that city
    let shortestPaths = new ConcurrentDictionary<City, int * City list>()

    let rec searchForShortestPath curCity distanceSoFar citiesVisitedSoFar =

        // Visit all available cities from the current city
        let visitAvailableDestinations() =
            let availableDestinations = distTweenCities.[curCity]

            // Loop through destinations and spawn new tasks
            for dest in availableDestinations.Keys do
                Task.Factory.StartNew(
                    new Action(fun () ->
                        searchForShortestPath
                            dest
                            (distTweenCities.[curCity].[dest] + distanceSoFar)
                            (citiesVisitedSoFar @ [dest])
                )) |> ignore

        // Have I already found a way to travel to this city?
        let visitedBefore = shortestPaths.ContainsKey(curCity)

        if not visitedBefore then
            // First time visiting, add to visted cities
            shortestPaths.TryAdd(curCity, (distanceSoFar, citiesVisitedSoFar))

            lock stdout
                (fun () -> printfn
                                "Original route to %A (%d) via: %A"
                                curCity distanceSoFar citiesVisitedSoFar)

            visitAvailableDestinations()

        else // We have visited this city before, let's see if this route is faster
            let shortestKnownPath, cities = shortestPaths.[curCity]

            if distanceSoFar < shortestKnownPath then
                // Update shortest path, revisit neighboring cities
                shortestPaths.[curCity] <- (distanceSoFar, citiesVisitedSoFar)

                lock stdout
                    (fun () -> printfn
                                    "Found shorter route to %A (%d) via: %A"
                                    curCity distanceSoFar citiesVisitedSoFar)

            visitAvailableDestinations()
        

    // Create the master task to find the shortest path between the two cities
    let t =
        Task.Factory.StartNew(
            new Action(fun () -> searchForShortestPath startingCity 0 [])
        )
    t.Wait()

    let dist, path = shortestPaths.[finalDestination]

    printfn
        "The shortest distance from %A to %A is %d miles, with route:\n%A"
        startingCity finalDestination
        dist path