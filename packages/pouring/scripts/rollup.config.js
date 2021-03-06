import { nodeResolve } from '@rollup/plugin-node-resolve';
import commonjs from '@rollup/plugin-commonjs';
import sourceMaps from 'rollup-plugin-sourcemaps';
import typescript from 'rollup-plugin-typescript2';



export default {
    input: `source/index.ts`,
    output: [
        {
            file: './distribution/index.js',
            format: 'cjs',
            sourcemap: true,
            exports: 'named',
        },
        {
            file: './distribution/index.es.js',
            format: 'es',
            sourcemap: true,
            exports: 'named',
        },
    ],
    external: [
        'commander',
        'cross-fetch',
        'encoding',
        'sync-fetch',
    ],
    watch: {
        include: 'source/**',
    },
    plugins: [
        nodeResolve(),
        commonjs(),
        sourceMaps(),
        typescript({
            file: '../tsconfig.json',
            useTsconfigDeclarationDir: true,
        }),
    ],
};
