require.config({
    paths: {
        velocity: "https://cdn.jsdelivr.net/velocity/1.2.3/velocity.min",
        interact: "https://cdn.jsdelivr.net/interact.js/1.2.6/interact.min",
        angular: "https://cdnjs.cloudflare.com/ajax/libs/angular.js/1.5.8/angular.min"
       },
    shim: {
        'angular': {
            exports: 'angular'
        }
    }
});
